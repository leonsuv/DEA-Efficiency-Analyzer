from cx_Freeze import setup, Executable
import sys

# Append build_exe command if no command is given.
if len(sys.argv) == 1:
    sys.argv.append("build_exe")

# Set base to "Win32GUI" on Windows to hide the console window.
base = None
if sys.platform == "win32":
    base = "Win32GUI"  # Change to None if you need a console

# Options for the build_exe command
build_exe_options = {
    "packages": [
         "os", "sys", "numpy", "pandas", "scipy", "nicegui", 
         "matplotlib", "dealib", "uvicorn", "starlette"
    ],
    "include_files": []  # Add any additional files you need here
}

setup(
    name="DEA Efficiency Analyzer",
    version="1.0",
    description="Executable for DEA Efficiency Analyzer",
    options={"build_exe": build_exe_options},
    executables=[Executable("main.py", base=base)]
)