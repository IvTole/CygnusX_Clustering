from pathlib import Path


def cube_data_path() -> Path:
    """
    Returns the location of the data cube, allowing for script executions in subfolders without worrying about the
    relative location of the data

    :return: the path to the data cube
    """
    cwd = Path("..")
    for folder in (cwd, cwd / "..", cwd / ".." / ".."):
        data_folder = folder / "data"
        if data_folder.exists() and data_folder.is_dir():
            print("Data (main) directory found in ", data_folder)
            return data_folder
        else:
            raise Exception("Data not found")
        
def plot_data_path() -> Path:
    """
    Returns the location of the plot directory, allowing for script executions in subfolders without worrying about the
    relative location of the data

    :return: the path to the plot directory
    """
    cwd = Path("..")
    for folder in (cwd, cwd / "..", cwd / ".." / ".."):
        data_folder = folder / "plots"
        if data_folder.exists() and data_folder.is_dir():
            print("Plot directory found in ", data_folder)
            return data_folder
        else:
            raise Exception("Plots directory not found")
        
def fits_data_path() -> Path:
    """
    Returns the location of the fits directory, allowing for script executions in subfolders without worrying about the
    relative location of the data

    :return: the path to the plot directory
    """
    cwd = Path("..")
    for folder in (cwd, cwd / "..", cwd / ".." / ".."):
        data_folder = folder / "fitsfiles"
        if data_folder.exists() and data_folder.is_dir():
            print("Fits directory found in ", data_folder)
            return data_folder
        else:
            raise Exception("Fits files directory not found")