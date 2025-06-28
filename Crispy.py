import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QTextEdit, QLineEdit, QLabel, QFileDialog
from PyQt5.QtCore import Qt, QRect
from PyQt5.QtGui import QFont, QColor, QPalette, QLinearGradient
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import re
import random

# Initialize codon table for standard genetic code
codon_table = CodonTable.unambiguous_dna_by_name["Standard"]

# Helper function to find synonymous codons
def get_synonymous_codons(codon):
    """Return list of synonymous codons for a given codon."""
    amino_acid = str(Seq(codon).translate())
    return [c for c in codon_table.forward_table if codon_table.forward_table[c] == amino_acid and c != codon]

# Main application window
class CRISPRWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("GeneSync : CRISPR Evader <beta version>")
        self.setFixedSize(800, 500)
        self.setWindowFlags(Qt.FramelessWindowHint)  # Frameless window
        self.init_ui()

    def init_ui(self):
        """Set up the modern, biotech-inspired GUI with a sleek and vibrant design."""
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout()
        layout.setSpacing(15)
        layout.setContentsMargins(20, 20, 20, 20)
        main_widget.setLayout(layout)

        palette = QPalette()
        gradient = QLinearGradient(0, 0, 0, 500)
        gradient.setColorAt(0, QColor(240, 245, 250))
        gradient.setColorAt(1, QColor(220, 230, 240))
        palette.setBrush(QPalette.Window, gradient)
        palette.setColor(QPalette.Button, QColor(0, 150, 136))
        palette.setColor(QPalette.ButtonText, QColor(255, 255, 255))
        self.setPalette(palette)
        main_widget.setStyleSheet("background: transparent;")

        title_bar = QWidget()
        title_bar.setFixedHeight(40)
        title_bar.setStyleSheet("""
            background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #26A69A, stop:1 #4FC3F7);
            border-radius: 5px;
        """)
        title_layout = QHBoxLayout()
        title_layout.setContentsMargins(10, 0, 10, 0)
        title_label = QLabel("GeneSync : Crispy <beta V-1.02>")
        title_label.setFont(QFont("Helvetica", 12, QFont.Bold))
        title_label.setStyleSheet("color: white;")
        close_btn = QPushButton("✕")
        close_btn.setFixedSize(30, 30)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #EF5350;
                color: white;
                border-radius: 15px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #F44336;
            }
        """)
        close_btn.clicked.connect(self.close)
        title_layout.addWidget(title_label)
        title_layout.addStretch()
        title_layout.addWidget(close_btn)
        title_bar.setLayout(title_layout)
        layout.addWidget(title_bar)

        seq_label = QLabel("Upload DNA Sequence (FASTA or Text)")
        seq_label.setFont(QFont("Helvetica", 10, QFont.Bold))
        seq_label.setStyleSheet("color: #37474F;")
        layout.addWidget(seq_label)
        
        self.seq_input = QTextEdit()
        self.seq_input.setFixedHeight(100)
        self.seq_input.setStyleSheet("""
            QTextEdit {
                border: 1px solid #CFD8DC;
                border-radius: 8px;
                padding: 5px;
                background-color: white;
                font-family: Helvetica;
            }
            QTextEdit:focus {
                border: 1px solid #26A69A;
            }
        """)
        layout.addWidget(self.seq_input)

        load_btn = QPushButton("Load Sequence File")
        load_btn.setFixedHeight(40)
        load_btn.setStyleSheet("""
            QPushButton {
                background-color: #26A69A;
                color: white;
                border-radius: 8px;
                font-family: Helvetica;
                font-size: 12px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #00897B;
            }
        """)
        load_btn.clicked.connect(lambda: self.load_sequence())
        layout.addWidget(load_btn)

        guide_label = QLabel("Enter CRISPR Guide RNA (20 bp)")
        guide_label.setFont(QFont("Helvetica", 10, QFont.Bold))
        guide_label.setStyleSheet("color: #37474F;")
        layout.addWidget(guide_label)
        
        self.guide_input = QLineEdit()
        self.guide_input.setFixedHeight(40)
        self.guide_input.setStyleSheet("""
            QLineEdit {
                border: 1px solid #CFD8DC;
                border-radius: 8px;
                padding: 5px;
                background-color: white;
                font-family: Helvetica;
            }
            QLineEdit:focus {
                border: 1px solid #26A69A;
            }
        """)
        layout.addWidget(self.guide_input)

        mutate_btn = QPushButton("Mutate & Save")
        mutate_btn.setFixedHeight(40)
        mutate_btn.setStyleSheet("""
            QPushButton {
                background-color: #26A69A;
                color: white;
                border-radius: 8px;
                font-family: Helvetica;
                font-size: 12px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #00897B;
            }
        """)
        mutate_btn.clicked.connect(lambda: self.process_sequence())
        layout.addWidget(mutate_btn)

        output_label = QLabel("Modified Sequence")
        output_label.setFont(QFont("Helvetica", 10, QFont.Bold))
        output_label.setStyleSheet("color: #37474F;")
        layout.addWidget(output_label)
        
        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        self.output_text.setStyleSheet("""
            QTextEdit {
                border: 1px solid #CFD8DC;
                border-radius: 8px;
                padding: 5px;
                background-color: white;
                font-family: Helvetica;
            }
        """)
        layout.addWidget(self.output_text)

        self.status_label = QLabel("Status: Ready")
        self.status_label.setFont(QFont("Helvetica", 10))
        self.status_label.setStyleSheet("color: #455A64;")
        layout.addWidget(self.status_label)

    def load_sequence(self):
        """Load DNA sequence from a FASTA or text file."""
        file_name, _ = QFileDialog.getOpenFileName(self, "Open Sequence File", "", "FASTA Files (*.fasta *.fa);;Text Files (*.txt)")
        if file_name:
            try:
                with open(file_name, "r") as file:
                    try:
                        seq_record = next(SeqIO.parse(file, "fasta"))
                        self.seq_input.setText(str(seq_record.seq).upper())
                    except:
                        file.seek(0)
                        seq = file.read().strip().upper()
                        self.seq_input.setText(seq)
                self.status_label.setText("Status: Sequence loaded successfully")
            except Exception as e:
                self.status_label.setText(f"Status: Error loading file: {str(e)}")

    def process_sequence(self):
        """Detect CRISPR sites and apply silent mutations."""
        sequence = self.seq_input.toPlainText().strip().upper()
        guide = self.guide_input.text().strip().upper()
        
        if not sequence or not guide:
            self.status_label.setText("Status: Please provide both sequence and guide RNA")
            return
        if len(guide) != 20:
            self.status_label.setText("Status: Guide RNA must be 20 bp")
            return
        if not all(c in "ATCG" for c in sequence + guide):
            self.status_label.setText("Status: Invalid DNA sequence or guide RNA")
            return

        modified_seq, modified_count, _ = apply_silent_mutations(sequence, guide)
        
        html_output = ""
        for i, (orig, mod) in enumerate(zip(sequence, modified_seq)):
            if orig != mod:
                html_output += f'<span style="background-color: #B2DFDB; color: #00695C;">{mod}</span>'
            else:
                html_output += mod
        self.output_text.setHtml(html_output)
        
        self.status_label.setText(f"Status: {modified_count} CRISPR sites evaded")

        file_name, _ = QFileDialog.getSaveFileName(self, "Save Modified Sequence", "", "FASTA Files (*.fasta);;Text Files (*.txt)")
        if file_name:
            with open(file_name, "w") as file:
                file.write(">Modified_Sequence\n" + modified_seq)
            self.status_label.setText(f"Status: {modified_count} CRISPR sites evaded, saved to {file_name}")

def apply_silent_mutations(sequence, guide):
    """Find CRISPR sites and introduce silent mutations."""
    modified_seq = list(sequence)
    modified_count = 0
    debug_log = []
    
    for match in re.finditer(guide, sequence):
        start, end = match.start(), match.end()
        debug_log.append(f"Match found at {start}-{end}")
        
        codon_start = start - (start % 3) if start % 3 != 0 else start
        if codon_start < start:
            codon_start += 3
        debug_log.append(f"Adjusted codon start: {codon_start}")
        
        for i in range(codon_start, min(end, len(sequence)), 3):
            codon = sequence[i:i+3]
            if len(codon) != 3:
                debug_log.append(f"Skipping incomplete codon at {i}")
                continue
                
            syn_codons = get_synonymous_codons(codon)
            if not syn_codons:
                debug_log.append(f"No synonymous codons for {codon} at {i}")
                continue
                
            new_codon = random.choice(syn_codons)
            if str(Seq(codon).translate()) == str(Seq(new_codon).translate()):
                modified_seq[i:i+3] = list(new_codon)
                modified_count += 1
                debug_log.append(f"Mutated {codon} to {new_codon} at {i}")
            else:
                debug_log.append(f"Translation mismatch for {codon} to {new_codon} at {i}")
    
    return "".join(modified_seq), modified_count, debug_log

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = CRISPRWindow()
    window.show()
    sys.exit(app.exec_())