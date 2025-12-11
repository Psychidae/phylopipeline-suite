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
        
        # Find the binary dynamically (folder name might vary)
        found_bin = None
        for root, dirs, files in os.walk(TOOLS_DIR):
            if "iqtree2" in files:
                found_bin = os.path.join(root, "iqtree2")
                break
        
        if found_bin and os.path.exists(found_bin):
            # Move to target location if it's not already there
            if os.path.abspath(found_bin) != os.path.abspath(iqtree_path):
                shutil.move(found_bin, iqtree_path)
            
            subprocess.run(["chmod", "+x", iqtree_path], check=True)
            msg = f"IQ-TREE 2 installed to {iqtree_path}"
        else:
            # List files for debugging
            files_list = []
            for r, d, f in os.walk(TOOLS_DIR):
                for file in f:
                    files_list.append(os.path.join(r, file))
            return False, f"Error: Binary not found. Extracted contents: {files_list[:5]}..."
            
        # Cleanup (Remove subdirectories created by tar)
        if os.path.exists(tar_path):
            os.remove(tar_path)
            
        # Remove any other directories in TOOLS_DIR that are not the binary itself
        # (This is a bit risky, so let's just leave them or be very specific. 
        #  For now, just finding and moving the binary is enough stability.)
        
        return True, msg

    except Exception as e:
        return False, f"Failed to install IQ-TREE 2: {e}"
