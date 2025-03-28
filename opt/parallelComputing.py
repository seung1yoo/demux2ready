#!/usr/bin/python3
import sys
import subprocess
import getopt
import time
import threading
from typing import List, Dict, Optional, Tuple
from pathlib import Path
import logging
from datetime import datetime

# Command dictionary for directory operations
CMD_DIC: Dict[str, str] = {
	"cd": "Path.cwd().chdir('%s')",
	"mkdir": "Path('%s').mkdir(parents=True, exist_ok=True)",
}

class TaskLogger:
	def __init__(self, task_id: str, log_dir: Path):
		self.task_id = task_id
		self.log_dir = log_dir
		self.log_dir.mkdir(parents=True, exist_ok=True)
		
		# Create log files with timestamp
		timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
		self.stdout_log = self.log_dir / f"{task_id}_{timestamp}_stdout.log"
		self.stderr_log = self.log_dir / f"{task_id}_{timestamp}_stderr.log"
		
		# Setup logging
		self.logger = logging.getLogger(f"task_{task_id}")
		self.logger.setLevel(logging.INFO)
		
		# File handler for stdout
		stdout_handler = logging.FileHandler(self.stdout_log)
		stdout_handler.setLevel(logging.INFO)
		stdout_formatter = logging.Formatter('%(asctime)s - %(message)s')
		stdout_handler.setFormatter(stdout_formatter)
		self.logger.addHandler(stdout_handler)
		
		# File handler for stderr
		stderr_handler = logging.FileHandler(self.stderr_log)
		stderr_handler.setLevel(logging.ERROR)
		stderr_formatter = logging.Formatter('%(asctime)s - ERROR - %(message)s')
		stderr_handler.setFormatter(stderr_formatter)
		self.logger.addHandler(stderr_handler)
		
		# Console handler
		console_handler = logging.StreamHandler()
		console_handler.setLevel(logging.INFO)
		console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		console_handler.setFormatter(console_formatter)
		self.logger.addHandler(console_handler)

	def log_stdout(self, message: str) -> None:
		"""Log stdout message."""
		self.logger.info(message)

	def log_stderr(self, message: str) -> None:
		"""Log stderr message."""
		self.logger.error(message)

	def get_log_files(self) -> Tuple[Path, Path]:
		"""Get the paths of log files."""
		return self.stdout_log, self.stderr_log

class ParallelComputing:
	def __init__(self):
		self.idx: int = 0
		self.script_list: List[str] = []
		self.max_threads: Optional[int] = None
		self.engine: Optional[str] = None
		self.msg: str = ""
		self.log_dir = Path("task_logs")
		self.task_loggers: Dict[str, TaskLogger] = {}

	def init(self) -> None:
		"""Initialize the parallel computing environment."""
		options, args = getopt.getopt(sys.argv[1:], "t:")
		for op, p in options:
			if op == "-t":
				self.max_threads = int(p)
			else:
				print("Unknown option:", op)
		
		if not args:
			raise ValueError("Engine type (SA) must be specified")
			
		self.engine = args[0]
		self.script_list = sys.stdin.read().split("\n")

		if self.max_threads is None and self.engine == "SA":
			self.max_threads = self._get_cpu_count()

	def _get_cpu_count(self) -> int:
		"""Get the number of CPU cores available."""
		return int(subprocess.getoutput("grep -c processor /proc/cpuinfo").strip())

	def _get_line(self) -> Optional[str]:
		"""Get the next line from script list."""
		if self.idx >= len(self.script_list):
			return None
		
		line = self.script_list[self.idx].strip()
		print(self.script_list[self.idx])
		self.msg += self.script_list[self.idx] + "\n"
		self.idx += 1
		return line

	def _run_thread(self, script: str, task_id: str) -> None:
		"""Execute a single script in a thread with logging."""
		logger = self.task_loggers[task_id]
		
		try:
			# Run the command and capture output
			process = subprocess.Popen(
				script,
				shell=True,
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE,
				text=True,
				bufsize=1,
				universal_newlines=True
			)
			
			# Read stdout and stderr in real-time
			while True:
				output = process.stdout.readline()
				if output:
					logger.log_stdout(output.strip())
				
				error = process.stderr.readline()
				if error:
					logger.log_stderr(error.strip())
				
				# Check if process has finished
				if output == '' and error == '' and process.poll() is not None:
					break
			
			# Get return code
			return_code = process.poll()
			if return_code != 0:
				logger.log_stderr(f"Process exited with return code {return_code}")
				
		except Exception as e:
			logger.log_stderr(f"Error executing script: {str(e)}")
			raise
		
		time.sleep(0)

	def _run_standalone(self, script_list: List[str]) -> None:
		"""Execute scripts in parallel using threads."""
		if self.max_threads is None:
			raise ValueError("max_threads must be set")

		threads: List[threading.Thread] = []
		for i, script in enumerate(script_list):
			cmd = script.split()[0]
			cmd_arg = script.split()[-1]
			
			if cmd in CMD_DIC:
				eval(CMD_DIC[cmd] % cmd_arg)
				continue

			# Create unique task ID
			task_id = f"task_{len(self.task_loggers)}_{i}"
			self.task_loggers[task_id] = TaskLogger(task_id, self.log_dir)

			while True:
				active_threads = [t for t in threads if t.is_alive()]
				if len(active_threads) < self.max_threads:
					thread = threading.Thread(
						target=self._run_thread,
						args=(script, task_id)
					)
					thread.start()
					threads.append(thread)
					break
				time.sleep(0.1)

		for thread in threads:
			thread.join()

	def _process_scripts(self) -> None:
		"""Process scripts based on control commands."""
		while True:
			script_list: List[str] = []
			script = self._get_line()
			
			if not script:
				break
				
			if "#join" in script:
				break
			elif "#fork_" in script:
				engine = script.split()[0].split("_")[-1]
				self._select_fork(engine)
			elif "#begin_" in script:
				engine = script.split()[0].split("_")[-1]
				self._select_begin(engine)
			elif "#fork" in script:
				self._fork()
			elif "#begin" in script:
				self._begin()
			else:
				script_list.append(script)
				if self.engine == "SA":
					self._run_standalone(script_list)

	def _fork(self) -> None:
		"""Process scripts in fork mode."""
		script_list: List[str] = []
		while True:
			script = self._get_line()
			if not script or "#join" in script:
				break
			elif "#fork_" in script:
				engine = script.split()[0].split("_")[-1]
				self._select_fork(engine)
			elif "#begin_" in script:
				engine = script.split()[0].split("_")[-1]
				self._select_begin(engine)
			elif "#fork" in script:
				self._fork()
			elif "#begin" in script:
				self._begin()
			else:
				script_list.append(script)
		
		if self.engine == "SA":
			self._run_standalone(script_list)

	def _begin(self) -> None:
		"""Process scripts in begin mode."""
		while True:
			script_list: List[str] = []
			script = self._get_line()
			
			if not script or "#end" in script:
				break
			elif "#fork_" in script:
				engine = script.split()[0].split("_")[-1]
				self._select_fork(engine)
			elif "#begin_" in script:
				engine = script.split()[0].split("_")[-1]
				self._select_begin(engine)
			elif "#fork" in script:
				self._fork()
			elif "#begin" in script:
				self._begin()
			else:
				script_list.append(script)
				if self.engine == "SA":
					self._run_standalone(script_list)

	def _select_fork(self, engine: str) -> None:
		"""Process scripts with specific engine in fork mode."""
		script_list: List[str] = []
		while True:
			script = self._get_line()
			if not script or "#join" in script:
				break
			elif "#fork_" in script:
				engine = script.split()[0].split("_")[-1]
				self._select_fork(engine)
			elif "#begin_" in script:
				engine = script.split()[0].split("_")[-1]
				self._select_begin(engine)
			elif "#fork" in script:
				self._fork()
			elif "#begin" in script:
				self._begin()
			else:
				script_list.append(script)
		
		if engine == "SA":
			self._run_standalone(script_list)

	def _select_begin(self, engine: str) -> None:
		"""Process scripts with specific engine in begin mode."""
		while True:
			script_list: List[str] = []
			script = self._get_line()
			
			if not script or "#end" in script:
				break
			elif "#fork_" in script:
				engine = script.split()[0].split("_")[-1]
				self._select_fork(engine)
			elif "#begin_" in script:
				engine = script.split()[0].split("_")[-1]
				self._select_begin(engine)
			elif "#fork" in script:
				self._fork()
			elif "#begin" in script:
				self._begin()
			else:
				script_list.append(script)
				if engine == "SA":
					self._run_standalone(script_list)

	def run(self) -> None:
		"""Main execution method."""
		try:
			self.msg += f"----------------------\nStart : {time.strftime('%Y/%m/%d/%H:%M:%S')}\n----------------------\n"
			print(self.msg)
			
			self._process_scripts()
			
			self.msg += f"----------------------\nEnd : {time.strftime('%Y/%m/%d/%H:%M:%S')}\n----------------------\n"
			print(self.msg)
			
		except Exception as e:
			error_msg = f"Error\n-----------------------------------\n{str(e)}\n-----------------------------------\n"
			print(error_msg)
			self.msg += error_msg
			sys.exit(1)

def main():
	if len(sys.argv) == 1:
		print(f"Usage : python {sys.argv[0]} [-t nThrds] <SA> <-.standard input>")
		sys.exit(1)
	
	pc = ParallelComputing()
	pc.init()
	pc.run()

if __name__ == "__main__":
	main()
