import os
import sys
import platform
import subprocess
import shutil

def install_tools():
    """
    Check if tools are installed and install them if missing (Linux only).
    """
    # Only run on Linux (Streamlit Cloud)
    if platform.system() != "Linux":
        print(f"Skipping tool installation on {platform.system()}")
        return

    # Define tool paths
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tools_dir = os.path.join(base_dir, "tools")
    iqtree_path = os.path.join(tools_dir, "iqtree2")

    # Create tools directory if it doesn't exist
    if not os.path.exists(tools_dir):
        os.makedirs(tools_dir)

    # Check for IQ-TREE 2
    if not os.path.exists(iqtree_path):
        print("Installing IQ-TREE 2 for Linux...")
        try:
            # Download URL for IQ-TREE 2.4.0 Linux (Intel 64-bit)
            url = "https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz"
            tar_path = os.path.join(tools_dir, "iqtree.tar.gz")
            
            # Download
            subprocess.run(["curl", "-L", "-o", tar_path, url], check=True)
            
            # Extract
            subprocess.run(["tar", "-xzf", tar_path, "-C", tools_dir], check=True)
            
            # Move binary (extracted folder name varies, usually iqtree-2.4.0-Linux)
            extracted_dir = os.path.join(tools_dir, "iqtree-2.4.0-Linux")
            bin_path = os.path.join(extracted_dir, "bin", "iqtree2")
            
            if os.path.exists(bin_path):
                shutil.move(bin_path, iqtree_path)
                subprocess.run(["chmod", "+x", iqtree_path], check=True)
                print(f"IQ-TREE 2 installed to {iqtree_path}")
            else:
                print("Error: Could not find iqtree2 binary in extracted files.")
                
            # Cleanup
            if os.path.exists(tar_path):
                os.remove(tar_path)
            if os.path.exists(extracted_dir):
                shutil.rmtree(extracted_dir)

        except Exception as e:
            print(f"Failed to install IQ-TREE 2: {e}")
    else:
        print("IQ-TREE 2 is already installed.")
