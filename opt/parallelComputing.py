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
		
		# Remove any existing handlers
		self.logger.handlers = []
		
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
		
		# Prevent propagation to root logger
		self.logger.propagate = False
		
		# Task execution info
		self.start_time = None
		self.end_time = None
		self.return_code = None
		self.command = None

	def log_stdout(self, message: str) -> None:
		"""Log stdout message."""
		self.logger.info(message)

	def log_stderr(self, message: str) -> None:
		"""Log stderr message."""
		self.logger.error(message)

	def get_log_files(self) -> Tuple[Path, Path]:
		"""Get the paths of log files."""
		return self.stdout_log, self.stderr_log

	def get_execution_info(self) -> Dict:
		"""Get task execution information."""
		duration = None
		if self.start_time and self.end_time:
			duration = self.end_time - self.start_time
			
		return {
			'task_id': self.task_id,
			'command': self.command,
			'start_time': self.start_time,
			'end_time': self.end_time,
			'duration': duration,
			'return_code': self.return_code,
			'status': 'SUCCESS' if self.return_code == 0 else 'FAILED',
			'log_files': {
				'stdout': str(self.stdout_log),
				'stderr': str(self.stderr_log)
			}
		}

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
		logger.command = script
		logger.start_time = datetime.now()
		
		try:
			# Log the command being executed
			logger.log_stdout("=" * 80)
			logger.log_stdout(f"Executing command: {script}")
			logger.log_stdout("=" * 80)
			
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
			logger.return_code = process.poll()
			if logger.return_code != 0:
				logger.log_stderr(f"Process exited with return code {logger.return_code}")
				
			# Log command completion
			logger.log_stdout("-" * 80)
			logger.log_stdout(f"Command completed with return code: {logger.return_code}")
			logger.log_stdout("-" * 80)
				
		except Exception as e:
			logger.log_stderr(f"Error executing script: {str(e)}")
			logger.return_code = 1
			raise
		finally:
			logger.end_time = datetime.now()
		
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

	def _print_summary(self) -> None:
		"""Print summary of all task executions."""
		print("\n" + "=" * 80)
		print("Task Execution Summary")
		print("=" * 80)
		
		total_tasks = len(self.task_loggers)
		successful_tasks = sum(1 for logger in self.task_loggers.values() if logger.return_code == 0)
		failed_tasks = total_tasks - successful_tasks
		
		print(f"\nTotal Tasks: {total_tasks}")
		print(f"Successful: {successful_tasks}")
		print(f"Failed: {failed_tasks}")
		print(f"Success Rate: {(successful_tasks/total_tasks)*100:.1f}%")
		
		print("\nDetailed Task Information:")
		print("-" * 80)
		for task_id, logger in self.task_loggers.items():
			info = logger.get_execution_info()
			duration_str = f"{info['duration'].total_seconds():.2f}s" if info['duration'] else "N/A"
			print(f"\nTask ID: {info['task_id']}")
			print(f"Command: {info['command']}")
			print(f"Status: {info['status']}")
			print(f"Return Code: {info['return_code']}")
			print(f"Duration: {duration_str}")
			print(f"Log Files:")
			print(f"  - stdout: {info['log_files']['stdout']}")
			print(f"  - stderr: {info['log_files']['stderr']}")
		
		print("\n" + "=" * 80)

	def run(self) -> None:
		"""Main execution method."""
		try:
			self.msg += f"----------------------\nStart : {time.strftime('%Y/%m/%d/%H:%M:%S')}\n----------------------\n"
			print(self.msg)
			
			self._process_scripts()
			
			self.msg += f"----------------------\nEnd : {time.strftime('%Y/%m/%d/%H:%M:%S')}\n----------------------\n"
			print(self.msg)
			
			# Print execution summary
			self._print_summary()
			
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
