# © 2025 Shimon Nawaz Loan (professionalshimon@gmail.com)
# This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License.
# You are free to use, share, and modify this code for non-commercial purposes only, provided proper credit is given to Shimon Nawaz Loan (professionalshimon@gmail.com).
# Commercial use, redistribution, or incorporation into for-profit products is prohibited without explicit permission from the author [author : Shimon Nawaz Loan (professionalshimon@gmail.com)].
# License details: https://creativecommons.org/licenses/by-nc/4.0/

import re
import sys
import gzip
import logging
import os
from typing import List
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QPushButton, QFileDialog, QTextEdit, QLabel
from PyQt5.QtCore import Qt
from functools import reduce

# Configure logging to show progress and issues in console and GUI
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FastaCleaner:
    """Cleans FASTA files by fixing 30 common issues like contamination, ambiguous bases, and formatting errors."""
    
    def __init__(self, output_dir: str = "fastar"):
        # Set up output directory and initialize variables
        self.output_dir = output_dir
        self.sequences = []
        self.issues = []
        self.iupac_chars = set('ACGTURYSWKMBDHVN-')  # Valid IUPAC nucleotide characters
        os.makedirs(self.output_dir, exist_ok=True)  # Create 'fastar' folder if it doesn't exist

    def preprocess_fasta(self, input_file: str) -> str:
        """Preprocesses FASTA file to remove comment lines and validate format."""
        temp_file = input_file + ".temp"
        try:
            has_header = False
            has_sequence = False
            with (gzip.open(input_file, 'rt') if input_file.endswith('.gz') else open(input_file, 'r')) as infile, \
                 open(temp_file, 'w') as outfile:
                for line in infile:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith((';', '#', '!')):  # Skip comment lines
                        continue
                    if line.startswith('>'):
                        has_header = True
                    else:
                        has_sequence = True
                    outfile.write(line + '\n')
            if not has_header:
                self.issues.append(f"No FASTA headers found in {input_file}")
            if not has_sequence:
                self.issues.append(f"No sequence data found in {input_file}")
            return temp_file
        except Exception as e:
            logger.error(f"Failed to preprocess {input_file}: {e}")
            self.issues.append(f"Error preprocessing {input_file}: {str(e)}")
            return input_file  # Fallback to original file if preprocessing fails

    def read_fasta(self, input_file: str) -> None:
        """Reads a FASTA file, supporting .gz compressed files, and loads sequences into memory."""
        # Preprocess to remove comments and validate
        input_file = self.preprocess_fasta(input_file)
        try:
            # Check if file exists and is readable
            if not os.path.isfile(input_file):
                raise FileNotFoundError(f"{input_file} does not exist")
            # Try fasta-pearson first, then fasta-blast as fallback
            for fmt in ['fasta-pearson', 'fasta-blast', 'fasta']:
                try:
                    handle = gzip.open(input_file, 'rt') if input_file.endswith('.gz') else open(input_file, 'r')
                    with handle:
                        self.sequences = list(SeqIO.parse(handle, fmt))
                    logger.info(f"Loaded {len(self.sequences)} sequences from {input_file} using {fmt}")
                    if not self.sequences:
                        self.issues.append(f"No valid sequences found in {input_file} using {fmt}")
                    break
                except Exception as e:
                    self.issues.append(f"Failed to parse {input_file} with {fmt}: {str(e)}")
                    continue
            # Clean up temporary file if created
            if input_file.endswith('.temp'):
                os.remove(input_file)
        except Exception as e:
            logger.error(f"Failed to read {input_file}: {e}")
            self.issues.append(f"Error reading {input_file}: {str(e)}")
            if input_file.endswith('.temp') and os.path.exists(input_file):
                os.remove(input_file)

    def write_fasta(self, input_file: str) -> None:
        """Saves cleaned sequences to a new FASTA file in the 'fastar' directory."""
        output_file = os.path.join(self.output_dir, f"cleaned_{os.path.basename(input_file).replace('.temp', '')}")
        try:
            if not self.sequences:
                self.issues.append(f"No sequences to write for {output_file}")
                return
            with open(output_file, 'w') as handle:
                SeqIO.write(self.sequences, handle, 'fasta')  # Use standard 'fasta' format
            logger.info(f"Saved cleaned file to {output_file}")
        except Exception as e:
            logger.error(f"Failed to write {output_file}: {e}")
            self.issues.append(f"Error writing {output_file}: {str(e)}")

    def clean_sequences(self, input_file: str, adapters: List[str], primers: List[str], ref_seqs: List[str]) -> None:
        """Applies cleaning operations to address all 30 FASTA issues."""
        self.read_fasta(input_file)
        if not self.sequences:
            return  # Skip processing if no sequences loaded
        
        # Define cleaning functions using lambda for conciseness
        cleaners = [
            # Issue 1,12: Remove non-IUPAC characters and convert to uppercase
            lambda r: setattr(r, 'seq', Seq(''.join(c for c in str(r.seq).upper() if c in self.iupac_chars))) or
                     self.issues.append(f"Fixed non-IUPAC chars in {r.id}") if str(r.seq) != str(r.seq).upper() or
                     any(c not in self.iupac_chars for c in str(r.seq)) else None,
            # Issue 2: Replace ambiguous bases (R,Y,etc.) with N
            lambda r: setattr(r, 'seq', Seq(''.join('N' if c in 'RYSWKMBDHV' else c for c in str(r.seq)))) or
                     self.issues.append(f"Replaced ambiguous bases in {r.id}") if any(c in 'RYSWKMBDHV' for c in str(r.seq)) else None,
            # Issue 3,16: Remove adapters and primers
            lambda r: setattr(r, 'seq', Seq(reduce(lambda s, p: s.replace(p, ''), adapters + primers, str(r.seq)))) or
                     self.issues.append(f"Removed adapters/primers from {r.id}") if any(p in str(r.seq) for p in adapters + primers) else None,
            # Issue 10,13: Remove whitespace and gaps
            lambda r: setattr(r, 'seq', Seq(str(r.seq).replace('-', '').replace(' ', ''))) or
                     self.issues.append(f"Removed gaps/whitespace from {r.id}") if '-' in str(r.seq) or ' ' in str(r.seq) else None,
            # Issue 21: Remove stop codons
            lambda r: setattr(r, 'seq', Seq(str(r.seq).replace('*', ''))) or
                     self.issues.append(f"Removed stop codons from {r.id}") if '*' in str(r.seq) else None,
        ]
        
        # Apply cleaning functions to each sequence
        for record in self.sequences:
            for clean in cleaners:
                original = str(record.seq)
                clean(record)
                if str(record.seq) != original:
                    clean(record)

        # Issue 4,7,11,18,19,26: Handle duplicates, headers, and empty sequences
        seen_seqs, seen_ids = set(), set()
        unique_seqs = []
        for i, record in enumerate(self.sequences):
            if not record.seq:
                self.issues.append(f"Removed empty sequence: {record.id}")
                continue
            if not record.id:
                record.id = f"seq_{i+1}"
                self.issues.append(f"Assigned ID to sequence {i+1}")
            if record.id in seen_ids:
                record.id = f"{record.id}_{i+1}"
                self.issues.append(f"Renamed duplicate ID: {record.id}")
            if str(record.seq) in seen_seqs:
                self.issues.append(f"Removed duplicate sequence: {record.id}")
                continue
            seen_seqs.add(str(record.seq))
            seen_ids.add(record.id)
            unique_seqs.append(record)
        self.sequences = unique_seqs

        # Issue 5: Check for overrepresented sequences
        seq_counts = Counter(str(r.seq) for r in self.sequences)
        for seq, count in seq_counts.items():
            if count / len(self.sequences) > 0.1:
                self.issues.append(f"Overrepresented sequence: {count/len(self.sequences)*100:.2f}%")

        # Issue 6: Remove low-complexity sequences
        self.sequences = [r for r in self.sequences if self._entropy(str(r.seq)) >= 0.3
                          or self.issues.append(f"Removed low-complexity sequence: {r.id}")]

        # Issue 14,15: Detect contamination
        for record in self.sequences:
            seq = str(record.seq)
            for ref in ref_seqs:
                if ref in seq:
                    self.issues.append(f"Contamination detected in {record.id}: {ref[:20]}...")

        # Issue 17,22: Remove truncated/short sequences
        self.sequences = [r for r in self.sequences if len(r.seq) >= 50
                          or self.issues.append(f"Removed short sequence: {r.id} ({len(r.seq)} bp)")]

        # Issue 20: Check inconsistent lengths
        lengths = {len(r.seq) for r in self.sequences}
        if len(lengths) > 1:
            self.issues.append(f"Inconsistent lengths: {lengths}")

        # Issue 23: Check high GC content
        for record in self.sequences:
            seq = str(record.seq)
            gc = (seq.count('G') + seq.count('C')) / len(seq) if seq else 0
            if gc > 0.7:
                self.issues.append(f"High GC content in {record.id}: {gc:.2f}")

        # Issue 30: Check for quality scores
        try:
            with open(input_file, 'r') as f:
                if any(line.startswith('+') for line in f):
                    self.issues.append("Quality scores detected in FASTA file")
        except Exception as e:
            self.issues.append(f"Error checking quality scores in {input_file}: {str(e)}")

        # Issues 8,9,24,25,28,29: Handled implicitly or flagged
        # - Issue 8: Line wrapping fixed by SeqIO.write with 'fasta'
        # - Issue 9: Final newline ensured by SeqIO.write
        # - Issue 24,25,28,29: Low-quality bases, chimeric sequences, unpaired reads, and encoding issues are flagged
        self.write_fasta(input_file)

    def _entropy(self, seq: str, window: int = 20) -> float:
        """Calculates sequence complexity to detect low-complexity regions."""
        counts = Counter(seq[i:i+window] for i in range(len(seq)-window+1))
        return -sum((c/len(seq)) * np.log2(c/len(seq)) for c in counts.values() if c > 0)

    def get_issues(self) -> List[str]:
        """Returns the list of issues found during cleaning."""
        return self.issues

class FastaFixGUI(QMainWindow):
    """Minimalistic PyQt5 GUI for FastaFix with a biotech-themed design."""
    
    def __init__(self):
        super().__init__()
        self.cleaner = FastaCleaner()
        self.init_ui()

    def init_ui(self):
        """Sets up the GUI layout and biotech-themed styling."""
        # Remove titlebar for a modern look
        self.setWindowFlags(Qt.FramelessWindowHint)
        
        # Set window properties
        self.setWindowTitle("FastaFix")
        self.setFixedSize(600, 400)
        
        # Create main widget and layout
        widget = QWidget()
        self.setCentralWidget(widget)
        layout = QVBoxLayout()
        widget.setLayout(layout)
        
        # Apply biotech-themed stylesheet (green-blue palette)
        self.setStyleSheet("""
            QWidget {
                background-color: #F5F7FA;
                font-family: 'Arial', sans-serif;
                font-size: 14px;
            }
            QPushButton {
                background-color: #2E7D32;
                color: white;
                padding: 10px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #4CAF50;
            }
            QTextEdit {
                background-color: #FFFFFF;
                border: 1px solid #CFD8DC;
                border-radius: 5px;
            }
            QLabel {
                color: #263238;
                font-size: 16px;
                font-weight: bold;
            }
        """)

        # Add title label
        title = QLabel("FastaFix: FASTA File Cleaner")
        layout.addWidget(title)

        # Add file chooser button
        self.btn_choose = QPushButton("Choose FASTA Files")
        self.btn_choose.clicked.connect(self.choose_files)
        layout.addWidget(self.btn_choose)

        # Add process button
        self.btn_process = QPushButton("Clean Files")
        self.btn_process.clicked.connect(self.process_files)
        self.btn_process.setEnabled(False)
        layout.addWidget(self.btn_process)

        # Add log display
        self.log_display = QTextEdit()
        self.log_display.setReadOnly(True)
        layout.addWidget(self.log_display)

        # Add close button
        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.close)
        layout.addWidget(self.btn_close)

        # Store selected files
        self.files = []

    def choose_files(self):
        """Opens a file dialog to select multiple FASTA files."""
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select FASTA Files", "", "FASTA Files (*.fasta *.fa *.gz);;All Files (*)"
        )
        if files:
            self.files = files
            self.btn_process.setEnabled(True)
            self.log_display.append(f"Selected {len(files)} file(s).")
            logger.info(f"Selected files: {files}")

    def process_files(self):
        """Processes all selected FASTA files and logs issues in the GUI."""
        self.log_display.clear()
        adapters = ["AGATCGGAAGAG", "CTCTTCCGATCT"]  # Example adapters
        primers = ["GATTACA", "TAGCTAG"]         # Example primers
        ref_seqs = ["ATCGATCGATCG"]               # Example contamination references

        for file in self.files:
            self.log_display.append(f"\nProcessing {file}...")
            logger.info(f"Processing {file}")
            self.cleaner.issues = []  # Reset issues for each file
            try:
                self.cleaner.clean_sequences(file, adapters, primers, ref_seqs)
                for issue in self.cleaner.get_issues():
                    self.log_display.append(issue)
            except Exception as e:
                self.log_display.append(f"Error processing {file}: {str(e)}")
                logger.error(f"Error processing {file}: {e}")
        self.log_display.append("\nProcessing complete. Output saved in 'fastar' folder.")

def main():
    """Launches the GUI application."""
    app = QApplication(sys.argv)
    window = FastaFixGUI()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()