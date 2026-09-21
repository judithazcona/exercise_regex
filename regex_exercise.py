import re

log_lines = [
    "2024-01-15 10:02:11 INFO Server started on port 8080",
    "2024-01-15 10:03:47 ERROR Failed to connect to database",
    "2024-01-16 08:15:00 WARNING Disk usage at 85%",
    "2024-01-16 14:22:39 ERROR Timeout while fetching https://example.com/api",
    "2024-01-17 09:00:05 INFO User admin logged in from 192.168.1.10",
    "2024-01-17 11:41:18 DEBUG Cache cleared successfully",
    "2024-01-18 03:12:56 ERROR Connection refused from 192.168.1.55",
    "2024-01-18 23:59:02 INFO Backup completed in 42s",
]

# Task 1
print ("\nTASK 1 \n")

#1
for line in log_lines:
    if re.match(r"2024-01-16", line):
        print(" ", line)

#2
for line in log_lines:
    if re.search(r"ERROR", line) or re.search(r"WARNING", line):
        print(" ", line)

#3
for line in log_lines:
    found = re.findall(r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3}", line)
    if found:
        print(" ", found)

#4
for line in log_lines:
    if re.search(r"\d+s$", line):
        print(" ", line)

#5
for line in log_lines:
    if re.search(r"https?://", line):
        print(" ", line)

#TASK 2
print ("\nTASK 2 \n")

def reverse_complement(sequence):
    complem = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    result = ""
    for base in sequence:
        result = result + complem[base]

    # invert the string!!
    return result[::-1]

class SequencingRead:
    def __init__(self, read_id, sequence):
        self.read_id = read_id
        self.sequence = sequence

# ^ start, $ end, .* anywhere,!! sequence has the patron at the start, and the reverse compliment at the end?
    def matches_mid_pair(self, forward_mid, reverse_mid): 
        patron = "^" + forward_mid + ".*" + reverse_complement(reverse_mid) + "$" #example: "^AT.*TC$"
        if re.fullmatch(patron, self.sequence):
            return True
        else:
            return False

#get the one in the mid
    def trim_mid_pair(self, forward_mid, reverse_mid):
        if self.matches_mid_pair(forward_mid, reverse_mid):
            inicio = len(forward_mid)
            final = len(self.sequence) - len(reverse_mid)
            return self.sequence[inicio:final]

    def describe(self):
            return "SequencingRead " + self.read_id + " (" + str(len(self.sequence)) + " bp)" # example: "SequencingRead demo_1 (46 bp)"


#demo
r1 = SequencingRead("demo_1", "AGCTTCGA" + "N" * 20 + reverse_complement("TGCAGGTC"))
print(r1.describe())
print(r1.matches_mid_pair("AGCTTCGA", "TGCAGGTC"))  # True
print(r1.matches_mid_pair("CGATCGAT", "GCTAGCTA"))  # False
print(r1.trim_mid_pair("AGCTTCGA", "TGCAGGTC"))     # 20 x "N"

# TASK 3

print ("\nTASK 3 \n")

# pip install biopython in the temrinal!!

import csv
import gzip
import os
import Bio.SeqIO

class Demultiplexer:
    def __init__(self, fasta_path, mid_table_path):
        self.reads = self.load_reads(fasta_path)
        self.mid_table = self.load_mid_table(mid_table_path)

        self.assigned = {}
        for label, forward_mid, reverse_mid in self.mid_table:
            self.assigned[label] = []

        self.unassigned = []

    def load_reads(self, fasta_path):
        read_list = []
        handle = gzip.open(fasta_path, "rt")
        for record in Bio.SeqIO.parse(handle, "fasta"):
            new_read = SequencingRead(record.id, str(record.seq))
            read_list.append(new_read)
        handle.close()
        return read_list

    def load_mid_table(self, mid_table_path):
        table = []
        handle = open(mid_table_path)
        reader = csv.DictReader(handle, delimiter=";")  # semicolon
        for row in reader:
            label = row["SampleID"] + "_" + row["Description"]
            forward_mid = row["FBarcodeSequence"]
            reverse_mid = row["RBarcodeSequence"]
            table.append((label, forward_mid, reverse_mid))
        handle.close()
        return table

    def assign_reads(self):
        for read in self.reads:
            was_assigned = False

            for label, forward_mid, reverse_mid in self.mid_table:

                # orientt 1
                if not was_assigned and read.matches_mid_pair(forward_mid, reverse_mid):
                    trimmed = read.trim_mid_pair(forward_mid, reverse_mid)
                    new_read = SequencingRead(read.read_id, trimmed)
                    self.assigned[label].append(new_read)
                    was_assigned = True

                # orientt 2
                if not was_assigned and read.matches_mid_pair(reverse_mid, forward_mid):
                    trimmed = read.trim_mid_pair(reverse_mid, forward_mid)
                    new_read = SequencingRead(read.read_id, trimmed)
                    self.assigned[label].append(new_read)
                    was_assigned = True

            if not was_assigned:
                self.unassigned.append(read)

    def report(self):
        result = ""
        for label, reads in self.assigned.items():
            result = result + label + "\t" + str(len(reads)) + "\n"
        result = result + "unassigned\t" + str(len(self.unassigned))
        return result

    def write_fasta(self, output_dir):
        os.makedirs(output_dir, exist_ok=True)

        for label, reads in self.assigned.items():
            if len(reads) == 0:
                continue

            path = output_dir + "/" + label + ".fasta"
            out = open(path, "w")
            for read in reads:
                out.write(">" + read.read_id + "\n")
                out.write(read.sequence + "\n")
            out.close()


#run
demux = Demultiplexer("fishes.fna.gz", "fishes_MIDs.csv")
demux.assign_reads()
print(demux.report())
demux.write_fasta("demux_output")