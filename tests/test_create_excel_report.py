import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import Mock, patch
from lib.create_excel_report import DataLoader, ReportGenerator

@pytest.fixture
def sample_md5sum_data():
    return pd.DataFrame({
        'md5': ['abc123', 'def456'],
        'file_path': ['/path/to/file1.fastq.gz', '/path/to/file2.fastq.gz']
    })

@pytest.fixture
def sample_s3upload_data():
    return pd.DataFrame({
        'File_name': ['file1.fastq.gz', 'file2.fastq.gz'],
        'File_size(Bytes)': ['1000', '2000'],
        'S3_url': ['s3://bucket/file1', 's3://bucket/file2'],
        'Presigned_url': ['http://presigned1', 'http://presigned2'],
        'Create_date': ['2024-01-01', '2024-01-01'],
        'Expiry_date': ['2024-02-01', '2024-02-01']
    })

@pytest.fixture
def sample_fastp_data():
    return pd.DataFrame({
        'Samples': ['sample1', 'sample2'],
        'ReadsCount': [1000000, 2000000],
        'BasesCount': [100000000, 200000000],
        'BasesCountGb': [0.1, 0.2],
        'GCRate': [50.0, 51.0],
        'Q20BaseRate': [98.0, 97.0],
        'Q30BaseRate': [95.0, 94.0],
        'DupRate': [5.0, 6.0],
        'InsertSize': [300, 301],
        'R1MeanLen': [150, 151],
        'R2MeanLen': [150, 151],
        'GoodReadRate': [95.0, 94.0]
    })

class TestDataLoader:
    def test_create_md5sum_info(self, sample_md5sum_data, tmp_path):
        # Create temporary md5sum file
        md5sum_file = tmp_path / "md5sum.txt"
        sample_md5sum_data.to_csv(md5sum_file, sep='\t', index=False)
        
        loader = DataLoader()
        result = loader.create_md5sum_info(md5sum_file)
        
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ['File Name', 'MD5 Checksum']
        assert len(result) == 2
        assert result['File Name'].iloc[0] == 'file1.fastq.gz'
        assert result['MD5 Checksum'].iloc[0] == 'abc123'

    def test_create_download_info(self, sample_s3upload_data, tmp_path):
        # Create temporary s3upload file
        s3upload_file = tmp_path / "s3upload.tsv"
        sample_s3upload_data.to_csv(s3upload_file, sep='\t', index=False)
        
        loader = DataLoader()
        result = loader.create_download_info(s3upload_file)
        
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ['File Name', 'File Size (Bytes)', 'S3 URL', 
                                      'Presigned URL', 'Create Date', 'Expiry Date']
        assert len(result) == 2
        assert result['File Name'].iloc[0] == 'file1.fastq.gz'
        assert result['File Size (Bytes)'].iloc[0] == 1000.0

    def test_create_sequencing_info(self, sample_fastp_data, tmp_path):
        # Create temporary fastp file
        fastp_file = tmp_path / "fastp.tsv"
        sample_fastp_data.to_csv(fastp_file, sep='\t', index=False)
        
        loader = DataLoader()
        result = loader.create_sequencing_info(fastp_file)
        
        assert isinstance(result, pd.DataFrame)
        assert 'Sample ID' in result.columns
        assert 'Total Reads' in result.columns
        assert len(result) == 2
        assert result['Sample ID'].iloc[0] == 'sample1'
        assert result['Total Reads'].iloc[0] == 1000000

class TestReportGenerator:
    @pytest.fixture
    def mock_args(self):
        args = Mock()
        args.fastp = "fastp.tsv"
        args.md5sum = "md5sum.txt"
        args.s3upload = "s3upload.tsv"
        args.output = "output.xlsx"
        return args

    def test_define_report_type_all_files(self, mock_args, tmp_path):
        # Create all required files
        (tmp_path / "fastp.tsv").touch()
        (tmp_path / "md5sum.txt").touch()
        (tmp_path / "s3upload.tsv").touch()
        
        with patch('pathlib.Path.exists', return_value=True):
            generator = ReportGenerator(mock_args)
            assert generator.report_type == 'all'

    def test_define_report_type_sequencing_only(self, mock_args):
        with patch('pathlib.Path.exists', side_effect=lambda x: x == "fastp.tsv"):
            generator = ReportGenerator(mock_args)
            assert generator.report_type == 'sequencing'

    @patch('pandas.ExcelWriter')
    def test_create_all_report(self, mock_excel_writer, mock_args, 
                             sample_fastp_data, sample_s3upload_data, sample_md5sum_data):
        with patch('pandas.read_csv') as mock_read_csv:
            mock_read_csv.side_effect = [sample_fastp_data, sample_s3upload_data, sample_md5sum_data]
            
            generator = ReportGenerator(mock_args)
            generator.create_all_report()
            
            # Verify Excel writer was called
            mock_excel_writer.assert_called_once()
            
            # Verify sheets were created
            mock_writer = mock_excel_writer.return_value.__enter__.return_value
            assert 'Sequencing Info' in mock_writer.sheets
            assert 'Download Info' in mock_writer.sheets 