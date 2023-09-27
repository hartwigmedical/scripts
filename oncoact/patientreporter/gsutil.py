import re

# Regex for gs file locations: gs:://<bucket-name>/<blob-name>
GS_PATH_REGEX = re.compile(r'^gs://([a-z0-9._-]+)((?:/[a-zA-Z0-9_.-]+)*)$')


def get_bucket_and_blob_from_gs_path(gs_path: str) -> (str, str):
    """
    Returns the bucket name and the blob name from a given google storage path.

    :param gs_path: the path in the following format: (gs://<bucket-name>/<blob-name>).
    :return: a tuple containing both the bucket name and blob name: (bucket, blob).
    """
    match = GS_PATH_REGEX.match(gs_path)
    if not match:
        raise ValueError(f"'{gs_path}' is not a valid path!")

    groups = match.groups()
    bucket = groups[0]
    blob = groups[1][1:] if groups[1] else None  # the [1:] is to remove the leading slash from the blob name

    return bucket, blob


def get_file_name_from_blob(blob_name: str) -> str:
    """
    Given a blob, it will treat it as a traditional file path and return the file name (the name after the last "/").

    :param blob_name: the full blob name, including (pseudo) subdirectories
    :return: the blobs file name
    """
    split = blob_name.split('/')
    return split[-1]
