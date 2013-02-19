from Bio import bgzf
from struct import unpack

def _encode_seq(seq):
    decoder = [
            "=", "A", "C", "M", "G", "R", "S", "V",
            "T", "W", "Y", "H", "K", "D", "B", "N",
            ]
    encoder = {}
    for item in decoder:
        encoder[item] = decoder.index(item)

    enc = 0
    for nt in seq:
        enc = enc << 4 | encoder[nt]
    return "%x" % enc


def _decode_seq(seq_str):
    decoder = [
            "=", "A", "C", "M", "G", "R", "S", "V",
            "T", "W", "Y", "H", "K", "D", "B", "N",
            ]
    decoded_seq = ""

    for char in seq_str:
        nibble1 = (ord(char) & 0xF0) >> 4
        nibble2 = ord(char) & 0x0F
        decoded_seq += decoder[nibble1]
        decoded_seq += decoder[nibble2]
    return decoded_seq


class Bam(object):
    val_lengths = {  # lengths of BAM types
            "A": 1,
            "Z": 1,
            "c": 1,
            "C": 1,
            "s": 2,
            "S": 2,
            "i": 4,
            "I": 4,
            }
    val_types = {  # Map between types of BAM and types of struct module
            "A": "s",
            "Z": "s",
            "c": "b",
            "C": "B",
            "s": "h",
            "S": "H",
            "i": "i",
            "I": "I",
            }



    def __init__(self, filename, keep_plaintext_targets=False):
        self.filename = filename
        self.keep_plaintext_targets = keep_plaintext_targets
        self.bgzf_file = bgzf.open(filename)
        magic = self.bgzf_file.read(4)
        if magic != "BAM\x01":
            raise ValueError("File does not start with the BAM magic code"
                    "BAM\\x01")

        # Reads the header
        header_length = unpack("<i", self.bgzf_file.read(4))[0]

        header = self.bgzf_file.read(header_length).strip("\x00")
        header_lines = header.split("\n")
        self.header = {}
        for line in header_lines:
            line = line.strip()
            line = line.strip("@")
            fields = line.split("\t")

            header_tag = fields[0]
            if not self.keep_plaintext_targets and header_tag == "SQ":
                next
            else:
                try:
                    self.header[header_tag].append(tuple(fields))
                except KeyError:
                    self.header[header_tag] = [tuple(fields), ]

        # Parse
        self.targets = []
        num_targets = unpack("<i", self.bgzf_file.read(4))[0]

        for iii in xrange(num_targets):
            name_len = unpack("<i", self.bgzf_file.read(4))[0]
            target_name = self.bgzf_file.read(name_len).strip("\x00")

            target_len = unpack("<i", self.bgzf_file.read(4))[0]
            self.targets.append((
                    target_name,
                    target_len,
                    ))
        print "targets:"
        print self.targets

    def next(self):
        try:
            aln_size_str = self.bgzf_file.read(4)
            if len(aln_size_str) != 4:
                raise StopIteration
        except:
            raise StopIteration
        aln_size = unpack("<i", aln_size_str)[0]
        print "aln_size is %i" % aln_size
        aln = self.bgzf_file.read(aln_size)
        print [c.encode("hex") for c in aln]

        target_id = unpack("<i", aln[:4])[0]
        print "TID: %i" % target_id
        aln = aln[4:]

        pos = unpack("<i", aln[:4])[0]
        print "Pos: %i" % pos
        aln = aln[4:]

        # This is little endian, so the last byte first, i.e. reverse the
        # fields
        len_read_name = unpack("<B", aln[:1])[0]
        print "LenReadName %i" % len_read_name
        aln = aln[1:]
        mapq = unpack("<B", aln[:1])[0]
        print "Mapq %i" % mapq
        aln = aln[1:]
        bin = unpack("<H", aln[:2])[0]
        print "Bin  %i" % bin
        aln = aln[2:]

        flag_nc = unpack("<I", aln[:4])[0]
        print "flag_nc %i %x" % (flag_nc, flag_nc)
        aln = aln[4:]
        flag = (flag_nc & 0xffff0000) >> 16
        print "flag %i" % flag
        n_cigar_op = (flag_nc & 0x0000ffff)
        print "NCigarOps %i" % n_cigar_op

        seq_len = unpack("<i", aln[:4])[0]
        print "seq_len %i" % seq_len
        aln = aln[4:]

        next_ref_id = unpack("<i", aln[:4])[0]
        print "next_ref_id %i" % next_ref_id
        aln = aln[4:]

        next_pos = unpack("<i", aln[:4])[0]
        print "next_pos %i" % next_pos
        aln = aln[4:]

        template_len = unpack("<i", aln[:4])[0]
        print "template_len %i" % template_len
        aln = aln[4:]

        read_name = aln[:len_read_name].strip("\x00")
        print "read_name %s" % read_name
        aln = aln[len_read_name:]

        cigar_str = aln[:n_cigar_op * 4]
        print "cigar_str %s" % cigar_str
        aln = aln[n_cigar_op * 4:]

        seq = aln[:((seq_len + 1) // 2)]
        print "seq %s" % seq
        print "decoded to %s" % _decode_seq(seq)
        aln = aln[((seq_len + 1) // 2):]

        qual = aln[:seq_len]
        print "qual %s" % qual
        aln = aln[seq_len:]
        print "len of aln str is %i" % len(aln)


        tags = {}
        while len(aln) > 4:
            tag = aln[:2]
            aln = aln[2:]

            val_type = aln[:1]
            aln = aln[1:]

            try:
                val_length = self.val_lengths[val_type]
                val_type = self.val_types[val_type]
            except KeyError:  # It's a BAM array type
                # the true value type is the next byte.
                val_type = aln[:1]
                aln = aln[1:]

                elements = unpack("<i", aln[:4])
                aln = aln[4:]

                val_length = self.val_lengths[val_type] * elements
                # type = "<12i" if the array is 12 integers, for struct module
                val_type = str(elements) + self.val_types[val_type]
            value = unpack("<%s" % val_type, aln[:val_length])
            aln = aln[val_length:]

            if len(value) == 1:  
                # no point in having the value in a tuple if it's a scalar
                value = value[0]

            tags[tag] = value
            print (tag, value)
