from stream import AlignmentStream

class AlignmentFilter(object):
    def __init__(self):
        self.edit_ratio = edit_ratio
        self.min_mismatch = min_mismatch
        self.coverage = coverage
        self.is_del = False
        self.is_unmapped = True
        self.is_duplicate = False

    def filter_with_stream(self):
        alignment = Alignment()
        for record in alignment:
            if self.coverage > record['coverage'] \
               and self.min_mismatch < record['mismatches']:
                pass
