import subprocess
import sys
import pkg_resources

def check_and_install_pip_packages(required_packages):
    for package in required_packages:
        try:
            # Check if package is installed; this will throw an exception if not
            pkg_resources.require(package)
        except pkg_resources.DistributionNotFound:
            print(f"{package} not found, installing...")
            install_with_pip(package)
        except pkg_resources.VersionConflict:
            pass  # Here, you could handle a version conflict or simply ignore it

def install_with_pip(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

def upgrade_with_pip(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", package])

def install_with_apt_get(package):
    subprocess.check_call(["sudo", "apt-get", "install", "-y", package])

# List of required pip packages, now including matplotlib and scikit-plot explicitly
required_pip_packages = [
    "graphviz",
    "matplotlib",  # Now explicitly listed for installation
    "numpy",
    "pandas",
    "joblib",
    "configparser",
    "scikit-learn",
    "scipy==1.11.4",
    "rpy2==3.5.15",
    "cffi",
    "scikit-plot"  # Added scikit-plot
]

# Upgrade pip before installing the packages
upgrade_with_pip("pip")

# Check and install required pip packages
check_and_install_pip_packages(required_pip_packages)

# Upgrade cffi and setuptools, and install wheel
upgrade_with_pip("cffi setuptools wheel")

# Install system-level dependency
install_with_apt_get("libffi-dev")

print("Installation of specified Python packages and system-level libraries is complete.")
