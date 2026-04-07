#!/usr/bin/env python3
import sys
import os

# Add current directory to path so onedgr package is found
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from onedgr.main import main_cli

if __name__ == "__main__":
    main_cli()
