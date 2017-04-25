#! /usr/bin/env python3

import vcf

__author__ = 'Shelley Brauneis'


class Assignment2:
    def __init__(self):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)
        self.filename='AmpliseqExome.20141120.NA24385.vcf'  #set file name

    def get_average_quality_of_son(self):
        self.file = vcf.Reader(open(self.filename, 'r'))
        sum=0
        count=0
        for record in self.file:         #for each record in the file, the record quality is collected in sum and the number of records in count
            sum += record.QUAL
            count = count + 1
        AvQualSon = sum/count    #the average quality for the son is the sum divided by the count
        return int(AvQualSon)   #to make the output look nicer, I chose to change the average to an int

    def get_total_number_of_variants_of_son(self):
        self.file = vcf.Reader(open(self.filename, 'r'))
        variants=0
        for record in self.file: #count the number of variants in file
            variants = variants + 1
        return variants

    def get_variant_caller_of_vcf(self):
        self.file = vcf.Reader(open(self.filename, 'r'))
        for line in self.file._header_lines:    #for every line in the header
            if line.startswith('##source'):     #if the line starts with info about the caller (##source)
                caller=line.split(" - ")        #split line at the " - "
                caller =caller[1].split("\"")   #the line is split using \", a list is then returned
        return caller[0]

    def get_human_reference_version(self):
        self.file = vcf.Reader(open(self.filename, 'r'))
        for line in self.file._header_lines:    #same as get variant caller of vcf, but this time we want reference instead of source
            if line.startswith('##reference=file://'):  #however, we have to specify if we want only the one result
                reference=line.split(".")
                reference=reference[0].split("/")
        return reference[7]

    def get_number_of_indels(self):
        self.file = vcf.Reader(open(self.filename, 'r'))
        indels=0
        for record in self.file:
            if record.is_indel:     #using record.is_indel (true/false) we can search the file for indels and add them using indels += 1
                indels += 1
        return indels

    def get_number_of_snvs(self):
        self.file = vcf.Reader(open(self.filename, 'r'))
        snv=0
        for record in self.file:    #same as with the indels, we can use record.is_snp (true/false) to count them (snp=snv for the purpose of this count)
            if record.is_snp:
                snv += 1
        return snv

    def get_number_of_heterozygous_variants(self):
        self.file = vcf.Reader(open(self.filename, 'r'))
        NrHet=0
        for record in self.file:        #we can use record.num_het to count the number of heterozygous variants, like when I used record.QUAL
            NrHet += record.num_het
        return NrHet

    def print_summary(self):
        print("Results: (this may take a while, 7 results total)")
        print("1. Average quality of the son:", self.get_average_quality_of_son())
        print("2. Variant number:", self.get_total_number_of_variants_of_son())
        print("3. Variant caller:", self.get_variant_caller_of_vcf())
        print("4. Human reference version:", self.get_human_reference_version())
        print("5. Number of indels:", self.get_number_of_indels())
        print("6. Number of snvs:", self.get_number_of_snvs())
        print("7. Number of heterozygous variants:", self.get_number_of_heterozygous_variants())



if __name__ == '__main__':
    print("Assignment 2", __author__)
    assignment1 = Assignment2()
    assignment1.print_summary()
