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
            print("Data directory found in ", data_folder)
            return data_folder
        else:
            raise Exception("Data not found")