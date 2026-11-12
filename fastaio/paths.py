"""Small helper for building output paths next to an input file."""

import os


class PATH:
    """Build output file and folder names from an input file path.

    Example
    -------
    >>> p = PATH("data/sample1.fasta")
    >>> p.add_to_name("distance", extention=".csv")
    'data/sample1_distance.csv'
    >>> p.go_up(1).add_to_dir("result").create_folder()   # creates ./result
    """

    def __init__(self, path):
        self.path = str(path)
        self.dir = os.path.dirname(self.path)
        self.name = os.path.splitext(os.path.basename(self.path))[0]
        self.extention = os.path.splitext(self.path)[1]

    def create_folder(self, folder_path=None):
        """Create ``folder_path`` (default: the current ``self.dir``)."""
        os.makedirs(folder_path or self.dir or ".", exist_ok=True)
        return self

    def add_to_name(self, add_name, dir=None, extention=None, after=True):
        """Return a path in ``dir`` (default ``self.dir``) named ``<name>_<add_name><ext>``
        (or ``<add_name>_<name><ext>`` when ``after=False``)."""
        directory = dir if dir else self.dir
        ext = extention if extention else self.extention
        stem = f"{self.name}_{add_name}" if after else f"{add_name}_{self.name}"
        return os.path.join(directory, stem + ext)

    def add_to_dir(self, add_name, dir=None):
        """Move ``self.dir`` one level down into ``add_name`` (chainable)."""
        base = self.dir if dir is None else dir
        self.dir = os.path.join(base, add_name)
        return self

    def go_up(self, how_much):
        """Move ``self.dir`` up ``how_much`` levels (chainable)."""
        for _ in range(how_much):
            self.dir = os.path.dirname(self.dir)
        return self

    def __str__(self):
        return self.path
