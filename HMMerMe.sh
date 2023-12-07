#!/bin/bash

#chmod -R g+w /ebio/ag-jekely/share/Nonbilaterian_neuropeptide_project/HMMer_search/input_HMM_profiles
#chmod -R g+w /ebio/ag-jekely/share/Nonbilaterian_neuropeptide_project/HMMer_search/input_Transcriptomes
umask 007


# level 1: go through the input transcriptome
for Transcriptome in input_Transcriptomes/*
do
	Filename="${Transcriptome#input_Transcriptomes/}"
	Species="${Filename%.fasta}"
	
	echo "starting with ${Species}"
	sed -i 's. ._ .g' "$Transcriptome"
	sed -i -e 's.\+._.g' -e 's.\:._.g' -e 's.\;._.g' -e 's.\/._.g' -e 's.|.___.g' -e 's.\\._.g' -e 's.{._.g' -e 's.\[._.g' -e 's.(._.g' -e 's.}._.g' -e 's.]._.g' -e 's.)._.g' -e 's.\^._.g' -e 's.\$._.g' -e 's.,._.g' -e 's.?._.g' -e 's.*..g' -e "s, ,\t,g" "$Transcriptome"

	mkdir tmp_conflicts_dir
	mkdir tmp_conflicts_dir/Identifiers_separated
	mkdir tmp_conflicts_dir/Identifiers_combined
	mkdir tmp_conflicts_dir/Sequences
	mkdir tmp_conflicts_dir/tmp_IDs
	mkdir tmp_conflicts_dir/Unique_identifiers
	mkdir tmp_conflicts_dir/Unique_tmp_IDs


	# level 2: search for the HMM profiles in the transcriptome
	for HMM_profile in input_HMM_profiles/*
	do 
		mkdir tmpoutput_dir

		HMMer_profile_filename="${HMM_profile#input_HMM_profiles/}"
		HMMer_geneclass="${HMMer_profile_filename%.hmm}"
		Prediction_file="${Species}__${HMMer_geneclass}_predictions.fasta"
		Identifier_list="${Species}__${HMMer_geneclass}_identifiers.fasta"
		Initial_HMMer_results="${Species}__${HMMer_geneclass}_evalues.fasta"
		Output1="${Species}__${HMMer_geneclass}_output1.fasta"
		Output2="${Species}__${HMMer_geneclass}_output2.fasta"
		Output3="${Species}__${HMMer_geneclass}_output3.fasta"
		Output4="${Species}__${HMMer_geneclass}_output4.fasta"
		
		echo "search for ${HMMer_geneclass}"

		# for setting different E-value, e.g 1e-5: add "-E 1e-5"
		hmmsearch -E 1e-10 "${HMM_profile}" "${Transcriptome}" > "tmpoutput_dir/${Output1}"
		
		sed -i '/--- full sequence ---/,/Domain annotation for each sequence (and alignments)/!d' "tmpoutput_dir/${Output1}"

		cp "tmpoutput_dir/${Output1}" "Nailed_it/original_HMMer_output/${Initial_HMMer_results}"	
		#sed -i 's/___/|/g' "Nailed_it/original_HMMer_output/${Initial_HMMer_results}"
		#sed -i 's/ //g' "Nailed_it/original_HMMer_output/${Initial_HMMer_results}"

		sed -i '4,/inclusion threshold/!d' "tmpoutput_dir/${Output1}"
		sed -i 's/ \+/	/g' "tmpoutput_dir/${Output1}"

		cut -f 10 "tmpoutput_dir/${Output1}" > "tmpoutput_dir/${Identifier_list}"
		
		cp "tmpoutput_dir/${Identifier_list}" "tmp_conflicts_dir/Identifiers_separated/${Identifier_list}"

		./.Apps/fastagrep.pl -f "tmpoutput_dir/${Identifier_list}" -X "${Transcriptome}" > "Nailed_it/original_HMMer_output/${Prediction_file}"
		cp "Nailed_it/original_HMMer_output/${Prediction_file}" "tmp_conflicts_dir/Sequences/${Prediction_file}"
		sed -i -e 's/___/|/g' -e "s,\t, ,g" "Nailed_it/original_HMMer_output/${Prediction_file}"
		sed -i 's._ . .g' "Nailed_it/original_HMMer_output/${Prediction_file}"

		rm -R -- tmpoutput_dir		
	done

	echo "--------------------------------------------------------------"
	echo " "
	echo "cleaning up"
	echo " "


	# check if the same sequences were retrieved with different HMM profiles (important when e.g. distinguishing between different GPCR classes)
	cat tmp_conflicts_dir/Identifiers_separated/* > tmp_conflicts_dir/Identifiers_combined/all_idents.fasta
	sort tmp_conflicts_dir/Identifiers_combined/all_idents.fasta > tmp_conflicts_dir/Identifiers_combined/all_idents_sorted.fasta
	uniq -d tmp_conflicts_dir/Identifiers_combined/all_idents_sorted.fasta > tmp_conflicts_dir/Identifiers_combined/all_double_idents.fasta
	uniq -u tmp_conflicts_dir/Identifiers_combined/all_idents_sorted.fasta > tmp_conflicts_dir/Unique_identifiers/uniq_identifiers.fasta

	# create separate ID files for IDs that were predicted once (no conflicts) and those that were predicted more than once (conflicts)
	for predicted_seqIDs in tmp_conflicts_dir/Identifiers_separated/*
	do
		predicted_seqIDs_file="${predicted_seqIDs#tmp_conflicts_dir/Identifiers_separated/}"
		predicted_seqIDs_type="${predicted_seqIDs_file%.fasta}"
		tmp1_conflict_seqsIDs="${predicted_seqIDs_type}_conflict_IDs_tmp1.fasta"
		tmp2_conflict_seqsIDs="${predicted_seqIDs_type}_conflict_IDs_tmp2.fasta"
		conflict_seqsIDs="${predicted_seqIDs_type}_conflict_IDs.fasta"
		tmp1_non_conflict_seq_IDs="${predicted_seqIDs_type}_non_conflict_IDs_tmp1.fasta"
		tmp2_non_conflict_seq_IDs="${predicted_seqIDs_type}_non_conflict_IDs_tmp2.fasta"
		non_conflict_seqIDs="${predicted_seqIDs_type}_no_conflict.fasta"


		cat "${predicted_seqIDs}" 'tmp_conflicts_dir/Identifiers_combined/all_double_idents.fasta' > "tmp_conflicts_dir/tmp_IDs/${tmp1_conflict_seqsIDs}"
		sort "tmp_conflicts_dir/tmp_IDs/${tmp1_conflict_seqsIDs}" > "tmp_conflicts_dir/tmp_IDs/${tmp2_conflict_seqsIDs}"
		uniq -d "tmp_conflicts_dir/tmp_IDs/${tmp2_conflict_seqsIDs}" > "tmp_conflicts_dir/tmp_IDs/${conflict_seqsIDs}"
		mv "tmp_conflicts_dir/tmp_IDs/${conflict_seqsIDs}" "Nailed_it/conflicted_HMMer_output/${conflict_seqsIDs}"

		cat "${predicted_seqIDs}" 'tmp_conflicts_dir/Unique_identifiers/uniq_identifiers.fasta' > "tmp_conflicts_dir/tmp_IDs/${tmp1_non_conflict_seq_IDs}"
		sort "tmp_conflicts_dir/tmp_IDs/${tmp1_non_conflict_seq_IDs}" > "tmp_conflicts_dir/tmp_IDs/${tmp2_non_conflict_seq_IDs}"
		uniq -d "tmp_conflicts_dir/tmp_IDs/${tmp2_non_conflict_seq_IDs}" > "tmp_conflicts_dir/tmp_IDs/${non_conflict_seqIDs}"
		mv "tmp_conflicts_dir/tmp_IDs/${non_conflict_seqIDs}" "Nailed_it/${non_conflict_seqIDs}"

		sed -i -e 's/___/|/g' -e "s,\t, ,g" "Nailed_it/conflicted_HMMer_output/${conflict_seqsIDs}"
		sed -i 's._ . .g' "Nailed_it/conflicted_HMMer_output/${conflict_seqsIDs}"
		sed -i -e 's/___/|/g' -e "s,\t, ,g" "Nailed_it/${non_conflict_seqIDs}"
		sed -i 's._ . .g' "Nailed_it/${non_conflict_seqIDs}"

	done
	
	rm -R -- tmp_conflicts_dir/tmp_IDs


	# get separate files with sequences that were only retrieved once (non conflicted) and those that were retrieved more thatn once (conflicted)
	for predicted_seqs in tmp_conflicts_dir/Sequences/*
	do
		predicted_seq_file="${predicted_seqs#tmp_conflicts_dir/Sequences/}"
		predicted_seq_type="${predicted_seq_file%.fasta}"
		conflict_seqs="${predicted_seq_type}_conflicts.fasta"
		non_conflict_seqIDs="${predicted_seq_type}_no_conflicts.fasta"

		./.Apps/fastagrep.pl -f tmp_conflicts_dir/Identifiers_combined/all_double_idents.fasta -X "${predicted_seqs}" > "Nailed_it//${conflict_seqs}"
		./.Apps/fastagrep.pl -f tmp_conflicts_dir/Unique_identifiers/uniq_identifiers.fasta -X "${predicted_seqs}" > "Nailed_it/${non_cconflicted_HMMer_outputonflict_seqIDs}"

		sed -i -e 's/___/|/g' -e "s,\t, ,g" "Nailed_it/conflicted_HMMer_output/${conflict_seqs}"
		sed -i 's._ . .g' "Nailed_it/conflicted_HMMer_output/${conflict_seqs}"
		sed -i -e 's/___/|/g' -e "s,\t, ,g" "Nailed_it/${non_conflict_seqIDs}"
		sed -i 's._ . .g' "Nailed_it/${non_conflict_seqIDs}"

	done


	rm -R -- tmp_conflicts_dir
	
	# transfer the processed transcriptome to the "processed folder"
	sed -i 's/___/|/g' "${Transcriptome}"
	sed -i 's._ . .g' "${Transcriptome}"
	mv -- "${Transcriptome}" "./Processed_input_files/Transcriptomes/${Filename}"



done

mv -- input_HMM_profiles/* Database_with_HMM_profiles/


grep -c -r ">" ./Nailed_it/original_HMMer_output/*_predictions.fasta > ./Nailed_it/original_HMMer_output/00_Number_of_all_predicted_genes.txt
sed -i "s,./Nailed_it/original_HMMer_output/,,g; s,__,\t,g; s,_predictions.fasta:,\t,g" ./Nailed_it/original_HMMer_output/00_Number_of_all_predicted_genes.txt

grep -c -r ">" ./Nailed_it/conflicted_HMMer_output/*_conflicts.fasta > ./Nailed_it/conflicted_HMMer_output/00_Number_of_conflicted_predicted_genes.txt
sed -i "s,./Nailed_it/conflicted_HMMer_output/,,g; s,__,\t,g; s,_conflicts.fasta:,\t,g" ./Nailed_it/conflicted_HMMer_output/00_Number_of_conflicted_predicted_genes.txt

grep -c -r ">" ./Nailed_it/*_no_conflicts.fasta > ./Nailed_it/00_Number_of_predicted_genes_without_conflicts.txt
sed -i "s,./Nailed_it/,,g; s,__,\t,g; s,_no_conflicts.fasta:,\t,g" ./Nailed_it/00_Number_of_predicted_genes_without_conflicts.txt



echo "--------------------------------------------------------------"
echo " "
echo "done!"
echo " "
echo ".hmm input profiles were moved to 'Database_with_HMM_profiles'"
echo " "
echo "input Transcriptomes were moved to 'Processed_input_files/Transcriptomes'"
echo " "
echo "original HMMer predictions are in 'Nailed_it/original_HMMer_output'"
echo " "
echo "sequences that were predicted by more than one HMM profile (= conflicts) are in 'Nailed_it/conflicted_HMMer_output'"
echo " "
echo "sequences that were only predicted once (= no conflicts) are in 'Nailed_it'"
echo " "
echo "number of predicted sequences (original, conflicts, no-conflicts) are listed in the 00_Number...txt files in the corresponding folders"
echo " "
echo "--------------------------------------------------------------"
