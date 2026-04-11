import os
import sys
import platform
import subprocess
import shutil
from modules.constants import TOOLS_DIR


def _is_valid_binary(path):
    """バイナリが実際に動作するかチェック"""
    if not os.path.exists(path):
        return False
    try:
        subprocess.run([path, "--version"], capture_output=True, timeout=10)
        return True
    except Exception:
        return False


def install_tools_linux():
    """
    Linux (Streamlit Cloud) 環境で IQ-TREE2 をダウンロード・インストールする。
    macOS では shutil.which でシステムインストール済みのツールを探す。
    戻り値: (success: bool, message: str)
    """
    system = platform.system()

    if system == "Darwin":
        # macOS: システムまたは tools/ にあれば OK
        for name in ["iqtree2", "iqtree"]:
            if shutil.which(name):
                return True, f"IQ-TREE found via PATH: {shutil.which(name)}"
        local = os.path.join(TOOLS_DIR, "iqtree2")
        if _is_valid_binary(local):
            return True, f"IQ-TREE found in tools/: {local}"
        return True, "macOS: IQ-TREE not found in PATH or tools/ — install manually or via Homebrew"

    if system != "Linux":
        return True, f"Skipping installer on {system}"

    # --- Linux インストール ---
    iqtree_path = os.path.join(TOOLS_DIR, "iqtree2")
    os.makedirs(TOOLS_DIR, exist_ok=True)

    if _is_valid_binary(iqtree_path):
        return True, "IQ-TREE 2 is ready."

    print("Installing IQ-TREE 2 for Linux...")

    # 古い無効バイナリを削除
    if os.path.exists(iqtree_path):
        try:
            os.remove(iqtree_path)
        except Exception:
            pass

    try:
        # アーキテクチャに応じたバイナリを選択
        machine = platform.machine()
        if machine in ("aarch64", "arm64"):
            url = "https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-arm.tar.gz"
        else:
            url = "https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz"

        tar_path = os.path.join(TOOLS_DIR, "iqtree.tar.gz")
        subprocess.run(["curl", "-fsSL", "-o", tar_path, url], check=True, timeout=120)
        subprocess.run(["tar", "-xzf", tar_path, "-C", TOOLS_DIR], check=True)

        # バイナリを探して移動
        found_bin = None
        for root, dirs, files in os.walk(TOOLS_DIR):
            if "iqtree2" in files:
                found_bin = os.path.join(root, "iqtree2")
                break

        if not found_bin:
            return False, f"iqtree2 binary not found after extraction in {TOOLS_DIR}"

        if os.path.abspath(found_bin) != os.path.abspath(iqtree_path):
            shutil.move(found_bin, iqtree_path)

        subprocess.run(["chmod", "+x", iqtree_path], check=True)

        # tarball とサブディレクトリを削除
        if os.path.exists(tar_path):
            os.remove(tar_path)

        return True, f"IQ-TREE 2 installed to {iqtree_path}"

    except subprocess.TimeoutExpired:
        return False, "IQ-TREE download timed out"
    except Exception as e:
        return False, f"Failed to install IQ-TREE 2: {e}"
