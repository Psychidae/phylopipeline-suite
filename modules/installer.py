import os
import sys
import platform
import subprocess
import shutil
from modules.constants import TOOLS_DIR

def install_tools_linux():
    """
    Check if tools are installed and install them if missing (Linux only).
    Returns (success, message)
    """
    # Only run on Linux (Streamlit Cloud)
    if platform.system() != "Linux":
        return True, f"Skipping installation on {platform.system()}"

    iqtree_path = os.path.join(TOOLS_DIR, "iqtree2")

    # Create tools directory if it doesn't exist
    if not os.path.exists(TOOLS_DIR):
        try:
            os.makedirs(TOOLS_DIR, exist_ok=True)
        except Exception as e:
            return False, f"Failed to create tools dir: {e}"

    # Function to check if binary is valid
    def is_valid_binary(path):
        if not os.path.exists(path):
            return False
        try:
            # Try running version command
            subprocess.run([path, "--version"], capture_output=True, check=True)
            return True
        except (OSError, subprocess.CalledProcessError):
            return False

    # Check for IQ-TREE 2
    if is_valid_binary(iqtree_path):
        return True, "IQ-TREE 2 is ready."

    print(f"Installing IQ-TREE 2 for Linux (Replacing invalid or missing binary)...")
    
    # Cleanup invalid binary if exists
    if os.path.exists(iqtree_path):
            try: os.remove(iqtree_path)
            except: pass
    
    try:
        # Download URL for IQ-TREE 2.4.0 Linux (Intel 64-bit)
        url = "https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz"
        tar_path = os.path.join(TOOLS_DIR, "iqtree.tar.gz")
        
        # Download
        subprocess.run(["curl", "-L", "-o", tar_path, url], check=True)
        
        # Extract
        subprocess.run(["tar", "-xzf", tar_path, "-C", TOOLS_DIR], check=True)
        
        # Move binary (extracted folder name varies, usually iqtree-2.4.0-Linux)
        extracted_dir = os.path.join(TOOLS_DIR, "iqtree-2.4.0-Linux")
        bin_path = os.path.join(extracted_dir, "bin", "iqtree2")
        
        if os.path.exists(bin_path):
            shutil.move(bin_path, iqtree_path)
            subprocess.run(["chmod", "+x", iqtree_path], check=True)
            msg = f"IQ-TREE 2 installed to {iqtree_path}"
        else:
            return False, "Error: Could not find iqtree2 binary in extracted files."
            
        # Cleanup
        if os.path.exists(tar_path):
            os.remove(tar_path)
        if os.path.exists(extracted_dir):
            shutil.rmtree(extracted_dir)
            
        return True, msg

    except Exception as e:
        return False, f"Failed to install IQ-TREE 2: {e}"
