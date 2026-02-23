import subprocess
import platform
#import importlib
import importlib.util

def check_nvidia_gpu():
    try:
        result = subprocess.run(["nvidia-smi"], capture_output=True, text=True)
        if result.returncode == 0:
            return True, result.stdout.split('\n')[2].strip()  # A line with GPU info
    except FileNotFoundError:
        pass
    return False, None

def check_apple_silicon_gpu():
    try:
        if platform.system() == "Darwin":
            result = subprocess.run(["system_profiler", "SPDisplaysDataType"], capture_output=True, text=True)
            if "Apple M" in result.stdout:
                return True, "Apple Silicon GPU (M1/M2/M3)"
    except Exception:
        pass
    return False, None

def check_gpu():
    print("Checking for available GPU...\n")

    nvidia, nvidia_info = check_nvidia_gpu()
    apple, apple_info = check_apple_silicon_gpu()

    if nvidia:
        print("NVIDIA GPU detected:")
        print(nvidia_info)
        return("nvidia")
    elif apple:
        print("Apple Silicon GPU detected:")
        print(apple_info)
        return("apple")
    else:
        print("No supported GPU detected (NVIDIA or Apple Silicon).")
        return("neither")

def check_module(module_name):
    spec = importlib.util.find_spec(module_name)
    return spec is not None

def check_cupy_or_mlx():
    print("Checking for GPU-related Python modules...\n")

    has_cupy = check_module("cupy")
    has_mlx = check_module("mlx")

    if has_cupy:
        print("'cupy' is installed (for NVIDIA GPUs).")
        return("cupy")
    elif has_mlx:
        print("'mlx' is installed (for Apple Silicon GPUs).")
        return("mlx")
    else:
        print("Neither 'cupy' or 'mlx' are installed")
        return("numpy")
