#!/usr/bin/env python3

import pandas as pd
import json
import numpy as np
import sys
import argparse


def errExit(message):
	sys.stderr.write("ERROR: %s\n"%message)
	sys.exit(1)

def main():

	#Map between domain names and fit_id in supplemental table
	domainToID= {"c-Src": 21879, "Fyn": 21873, "Grb2": 22075, "Lyn": 21882, "Yes": 21921, "Blk": 21915}

	#Builds the argument parser
	parser = argparse.ArgumentParser(description='Scores a protein sequence using a SH2 affinity model, outputting the predicted affinity sum.')
	domain_group = parser.add_mutually_exclusive_group(required=True)
	domain_group.add_argument('-d',   metavar='domain',     help='Name of domain %s'%list(domainToID.keys()), type=str, choices=domainToID.keys())
	domain_group.add_argument('-i',   metavar='fit_id',     help='Identity number of model, as listed in supplemental table', type=int)
	domain_group.add_argument('-j',   metavar='model.json', help='File containg a scoring model JSON object', type=str)
	sequence_group = parser.add_mutually_exclusive_group(required=True)
	sequence_group.add_argument('-s', metavar='sequence',   help='Specifies a protein sequence.')
	sequence_group.add_argument('-t', metavar="text_file",  help='Text file containing a list of protein sequences.')
	sequence_group.add_argument('-f', metavar="fasta_file", help='Fasta file containing the protein sequences.')
	parser.add_argument("--profile",                     help="Outputs the predictions for all offsets (instead of the sum)", action="store_true")
	args = parser.parse_args()

	#Reads the input sequences
	#OPTION: Provided sequence
	if args.s is not None:
		sequences = [[None, args.s]]
	#OPTION: Text file
	elif args.t is not None:
		sequences = []
		try:
			with open(args.t) as f:
				for l in f:
					sequences.append([None, l.rstrip()])	
		except FileNotFoundError:
			errExit("The input text file '%s' does not exist"%args.t)
	#OPTION: Fasta file
	elif args.f is not None:
		sequences = []
		try:
			with open(args.f) as f:
				for lName in f:
					seqName = lName.rstrip().lstrip(">")
					seq = f.readline().rstrip()
					sequences.append([seqName, seq])	
		except FileNotFoundError:
			errExit("The input text file '%s' does not exist"%args.t)
	else:
		errExit("Invalid input option.") #This can't happen

	#Reads the model
	domain = args.d
	if args.j is not None:
		#Reads model from file
		try:
			with open(args.j) as f:
				bmText = f.read()
		except FileNotFoundError:
			errExit("The JSON file '%s' does not exist"%args.f)

	else:
		#Reads the model from the supplement
		if args.d is not None:
			fit_id = domainToID[args.d]
		else:
			fit_id = args.i	

		df_models = pd.read_csv('S2_ProBound_models.tsv', sep='\t')
		#Selects model with correct fit_id
		filteredTable = df_models.loc[df_models["Fit ID"] == fit_id]
		#Checks so a model was located
		if len(filteredTable)==0:
			errExit("Cannot find model with fit_id=%d"%fit_id)
		#Gets the JSON in text format
		bmText = filteredTable.iloc[0]["Binding Mode 2 JSON"]

	bmJSON = json.loads(bmText)

	#Determines the alphabet and the order of the letters
	alphabet = bmJSON["modelSettings"]["letterOrder"]
	alphabetSet = set([c for c in alphabet])
	nAlpha = len(alphabet)
	charToI = dict([ (c,i) for i,c in  enumerate(alphabet)])

	#Gets the model
	size = bmJSON["modelSettings"]["bindingModes"][0]["size"]
	betas = bmJSON["coefficients"]["bindingModes"][0]["mononucleotide"]


	#Loops over sequences
	for (seqName, seq) in sequences:
		#Checks so the sequences are valid
		missingCharacters = set([c for c in seq if c not in alphabetSet])
		if len(missingCharacters)>0:
			errExit("The sequence '%s' contains the letter(s) '%s' that do not belong to the protein alphabet '%s'."%(seq, "".join(missingCharacters), alphabet))

		#Scores the sequence unless the sequence is shorter than the scoring matrix 
		if len(seq)>=size:
			#Scores the sequence at all offesets (x) by summing over all columns of the matrix (i)
			seqScores = [ np.exp(sum([betas[i*nAlpha+charToI[seq[x+i]]] for i in range(size)])) for x in range(len(seq)-size+1)]
		else:
			seqScores = None

		#Formats the output values
		if seqScores is None:
			output = "None"
		else:
			if args.profile:
				output = ",".join(["%e"%v for v in seqScores])
			else:
				output = "%e"%sum(seqScores)
	
		#Prints the output
		if seqName is None:
			print("%s\t%s"%(seq, output))
		else:
			print("%s\t%s"%(seqName, output))
		
main()
