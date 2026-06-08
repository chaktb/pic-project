"""
AWG WDM Design Utility - Main Application
Tkinter GUI port from Visual Basic 5.0/6.0

Original: ETRI, (C) 2000 by Jang-Uk Shin and ETRI
Python port: 2024
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import sys
import os

# Add project directory to path
sys.path.insert(0, os.path.dirname(__file__))

import globals as g
from frm_variables import FrmVariables

def main():
    root = tk.Tk()
    app = FrmVariables(root)
    root.mainloop()


if __name__ == "__main__":
    main()
