import requests, sys
import os
import pandas as pd
from Bio.Align import PairwiseAligner
import time
import shutil

def get_ensembl_transcripts_exons(ensembl):
    '''
    Fetch transcript and exon information from Ensembl API using the provided Ensembl gene ID.
    Handles retries and timeouts for connection issues.
    exon 0 in transcirpt 1
    data['Transcript'][1]['Exon'][0]['id'] - from ensmbl
    transcripts - list of all transcripts in the gene
    exons - for each transcript all its exons
    '''
    gene_data = {}
    base_url = f"https://rest.ensembl.org"
    #endpoint = f"/lookup/symbol/{species}/{gene_identifier}?expand=1" #use gene to get data
    endpoint = f"/lookup/id/{ensembl}?expand=1" #use ensembl id to get data
    headers = {"Content-Type": "application/json"}
    # Get and print Ensembl version
    # Retry logic parameters
    max_retries = 5
    retry_delay = 5  # Start with 5 seconds delay
    timeout = 10  # Timeout in seconds for each re
    version_endpoint = "/info/software"
    response = requests.get(base_url + version_endpoint, headers=headers, timeout=timeout)
    if response.ok:
        version_info = response.json()
        print("biomart version: ", version_info.get("release"))
    else:
        print("no version")
    for attempt in range(1, max_retries + 1):
        try:
    # Make a GET request to the Ensembl API
            response = requests.get(base_url + endpoint, headers=headers, timeout=timeout)
            # If the request was successful, exit the retry loop
            if response.status_code == 200:
                gene_data = response.json()
                break
            else:
                print(f"Attempt {attempt}: Received status code {response.status_code}. Retrying...")
        except requests.exceptions.Timeout:
            print(f"Attempt {attempt}: Request timed out. Retrying...")
        #Exponential backoff
        time.sleep(retry_delay)
        retry_delay *= 2  # Double the delay after each attempt
    
    # If we exhausted all retries and still failed    
    if response.status_code != 200: #i didn't get the response
    	return('','','','no_response','')
    #get canonical transcript data
    canonical_transcript = gene_data.get('canonical_transcript')
    canonical_exon_dict = {}
    if canonical_transcript:
        canonical_transcript_id = canonical_transcript.split('.')[0]
        canonical_exon_dict = get_canonical_transcript_exons(canonical_transcript_id)
    else:
        canonical_transcript_id = "no_canonical_transcript"
    # Extract transcript information from the response
    if "Transcript" in gene_data:
        transcript_info = []
        unique_exons = set()
        for transcript in gene_data.get("Transcript", []):
            for exon in transcript.get("Exon", []):
               exon_id = exon.get("id")
               if exon_id:
                   unique_exons.add(exon_id)
            info = {
                "id": transcript.get("id"),
                "is_canonical": transcript.get("is_canonical", False),
                "biotype": transcript.get("biotype"),
                "start": transcript.get("start"),
                "end": transcript.get("end"),
                "translation_start": transcript.get("Translation", {}).get("start"),
                "translation_end": transcript.get("Translation", {}).get("end"),
                "length": transcript.get("length"),
                "exons": transcript.get("Exon", [])
            }
            transcript_info.append(info)
        # Sort the transcripts
        sorted_transcripts = reorder_transcripts(transcript_info)
        chr = "chr" + gene_data['seq_region_name']
        if 'display_name' not in gene_data: #no gene name
            gene_data['display_name'] = "novel_gene" #in some ensmbl there are no gene names
        return (sorted_transcripts, gene_data['strand'], chr, gene_data['display_name'], canonical_exon_dict, unique_exons) 
    else: #no transcript
        return ([], gene_data['strand'], "chr", gene_data['display_name'], {}, []) 

def reorder_transcripts(transcripts):
    canonical = []
    protein_coding = []
    others = []
    for t in transcripts:
        if t.get("is_canonical") == 1:
            canonical.append(t)
        elif t.get("biotype") == "protein_coding":
            protein_coding.append(t)
        else:
            others.append(t)
    # sort by the existing 'length' field (descending)
    protein_coding_sorted = sorted(protein_coding, key=lambda x: x["length"], reverse=True)
    others_sorted = sorted(others, key=lambda x: x["length"], reverse=True)
    # final ordering
    return canonical + protein_coding_sorted + others_sorted

    
def get_canonical_transcript_exons(canonical_transcript_id, max_retries=10, delay=2):
    '''
    Parameters:
    canonical_transcript_id: Ensembl ID of the canonical transcript
    max_retries: Maximum number of retry attempts
    delay: Delay in seconds between retries
    Returns:
    exons_dict: a dict of the exons in the canonical transcript and their position (E1, E2, etc),
                or an error message if unsuccessful after retries.    
    '''
    # Get the exons of the canonical transcript
    url = f"https://rest.ensembl.org/lookup/id/{canonical_transcript_id}?expand=1&content-type=application/json"
    exon_dict = {}
    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                exons = data.get('Exon')
                if not exons:
                    print("No exons found for the canonical transcript")
                    return
                for index, exon in enumerate(exons):
                    if index == len(exons) - 1: #last exon is mark
                        exon_dict[exon['id']] = f"E{index+1}Last"
                    else:
                        exon_dict[exon['id']] = f"E{index+1}"
                #exon_dict = {exon['id']: f"E{index+1}" for index, exon in enumerate(exons)}
                return exon_dict
        except requests.exceptions.RequestException as e:
            print(f"Attempt {attempt+1}: Connection error: {e}. Retrying...")
        time.sleep(delay)  
    print("Error: Unable to fetch the exons after multiple attempts")
    return

def get_seq_exons(tra_exons_dict, biotype, strand):
    # type = genomic,cds,cdna,protein
    #get exon seq and start and end of each exon
    exons_duplication = {} #each exon one time
    exon_start_end = {}
    seq = {}
    server = "https://rest.ensembl.org"
    for index in range(0,len(tra_exons_dict)):
        #protein coding transcript
        start = tra_exons_dict[index]['translation_start']
        end = tra_exons_dict[index]['translation_end']
        if not isinstance(start, (int, float)): #not protein coding transcript
            start = tra_exons_dict[index]['start']
            end = tra_exons_dict[index]['end']
        exons_list = tra_exons_dict[index]['exons'] 
        for exon_index in range(0, len(exons_list)):
            name = exons_list[exon_index]['id']
            #print(tra_exons_dict[index]['id'],name)
            if name not in exons_duplication:
                #print(name)
                exons_duplication[name] = "y" #each exon get seq once
                ext = "/sequence/id/" + name + "?type=" + biotype 
                count_sec = 1
                max_retries = 15
                retries = 0
                while retries < max_retries:
                   try:
                       r = requests.get(server+ext, headers={"Content-Type": "text/x-fasta"}, timeout=10)
                       if r.ok: 
                           seq[name] = r.text
                           exon_start = int(r.text.split(":")[3])
                           exon_end = int(r.text.split(":")[4]) 
                           if exon_end > end or exon_start < start: #there is utr in the exon
                               seq[name] = get_cds(seq[name], start, end, exon_start, exon_end, strand)
                           exon_start_end[name] = r.text.split(":")[3] + "_" + r.text.split(":")[4]
                           break  # exit the loop if successful
                       else:
                           retries += 1
                           time.sleep(count_sec * 5)  # Exponential backoff between retries
                           count_sec += 1
                   except requests.exceptions.RequestException as e:
                        retries += 1
                        print(f"Attempt {retries} failed: {e}")
                        time.sleep(count_sec * 5)  # Exponential backoff between retries
                        count_sec += 1
                if retries == max_retries:
                    print(f"Failed to retrieve sequence for {name} after {max_retries} attempts.")
    return (seq, exon_start_end)

def get_cds(seq, start, end, exon_start, exon_end, strand):
    seq_name = seq.splitlines()[0]
    sequence = "".join(seq.splitlines()[1:])
    if exon_end > end and strand == 1:
        cds_seq = sequence[0:end-exon_start+1] 
    if exon_start < start and strand == 1:
        cds_seq = sequence[start-exon_start:]
    if exon_end > end and strand == -1:
        cds_seq = sequence[exon_end-end:]
    if exon_start < start and strand == -1:
        cds_seq = sequence[0:exon_end-start+1]
    if cds_seq:
       seq = seq_name + "\n" + cds_seq + "\n"
    return(seq) #if no cds in the exon, the original exon seq is return

def alignment(exons_seq_h, exons_seq_m, dir):
    concatenated_dict = {**exons_seq_h, **exons_seq_m}
    all_id = list(concatenated_dict)
    df_id = pd.DataFrame(index=all_id, columns=all_id)
    df_score = pd.DataFrame(index=all_id, columns=all_id)
    aligner = PairwiseAligner() 
    #copilot suggestion
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -12
    aligner.extend_gap_score = -4
    #aligner.target_end_gap_score = 0 - defualt - Enforces full-length alignment without free ends
    #aligner.query_end_gap_score = 0
    
    for x in range(0, len(concatenated_dict)):
        name_1 = all_id[x]
        name_1_seq = "\n".join(concatenated_dict[name_1].splitlines()[1:])#seq
        for y in range(0, len(concatenated_dict)):
            name_2 = all_id[y]
            name_2_seq = "\n".join(concatenated_dict[name_2].splitlines()[1:])#seq
            if len(name_1_seq) > 0 and len(name_2_seq) > 0:
                alignments = aligner.align(name_1_seq, name_2_seq)
                alignment = alignments[0]
                #print(alignment)
                identities = sum(aa1 == aa2 for aa1, aa2 in zip(alignment[0], alignment[1]))
                percentage_identity = round((identities / len(name_1_seq)) * 100)
                df_score.loc[name_1, name_2] = alignments.score
            else: #exons that are only UTR
                percentage_identity = -100
                df_score.loc[name_1, name_2] = -100
            df_id.loc[name_1, name_2] = percentage_identity

    output = dir + "/alignment_identity.txt"
    df_id.to_csv(output, sep="\t", index=True)
    output = dir + "/alignment_score.txt"
    df_score.to_csv(output, sep="\t", index=True)
    return(df_id)

def get_rank(canonical_exons, exon_id):
    rank = canonical_exons.get(exon_id, "*")
    return(rank)

def junction_position(position, strand, chr):
    #make sure junctions position is from the + strand
    if strand == 1:
        start = str(chr) + ':' + str(position[0])
        end = str(chr) + ':' + str(position[1])
    elif strand == -1:
        start = str(chr) + ':' + str(position[1])
        end = str(chr) + ':' + str(position[0])
    else:
        print('no_strand')
        return(0, 0)
    return(start, end)

def main():
  #/gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl113->qsub -cwd -V -q tals.q run_transcripts_exons_all_HN3.sh
    # Get the value from the command-line argument - the value is i, the gene index from the orthologs table
    version = 'HN6'
    min_id = 50 #min identity for orthologs exons
    # Check if a command-line argument is provided
    if len(sys.argv) != 2:
        print("Usage: python get_input.py <value>")
        sys.exit(1)
    i = int(sys.argv[1])
    working_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/'
    df_input = working_dir + 'all/good_ortologs_df.txt'
    filtered_df = pd.read_csv(df_input, sep='\t')
    columns_junctions = ['ensembl', 'symbol', 'exon_h', 'rank_h', 'position_h', 'exon_m', 'rank_m', 'position_m','%id_human','%id_mouse','Average %id']
    columns_statistics = ['ensembl','gene','transcripts','exons','exons_with_ortholog','E1','E2']

    ensembl_h = filtered_df['hsapiens_homolog_ensembl_gene'][i]
    ensembl_m = filtered_df['ensembl_gene_id'][i]
    gene = filtered_df['hsapiens_homolog_associated_gene_name'][i]
    #get the transcripts of the genes and for each transcript, its transclation start and end and its exons
    (tra_exons_h, strand_h, chr_h, gene_h, canonical_exons_h, all_exons_h) = get_ensembl_transcripts_exons(ensembl_h)
    if gene_h == 'no_response': 
        print(ensembl_h, 'no_response')
        sys.exit()
    if not isinstance(gene, float) and gene != gene_h:
        gene_h = gene_h + "," + gene #no gene name
    (tra_exons_m, strand_m, chr_m, gene_m, canonical_exons_m, all_exons_m) = get_ensembl_transcripts_exons(ensembl_m)
    if gene_m == 'no_response': 
        print(ensembl_h, ensembl_m, 'no_response')
        sys.exit()
    
    # Check number of exons
    num_exons_h = len(all_exons_h)
    num_exons_m = len(all_exons_m)
    if num_exons_h > 1 and num_exons_m > 1: #at least 2 exons
        exon_count_h = {}
        exon_count_h['E1'] = 0
        exon_count_h['E2'] = 0
        exon_count_m = {}
        exon_count_m['E1'] = 0
        exon_count_m['E2'] = 0
        print(i, gene_h, ensembl_h, ensembl_m)
        new_dir = working_dir + 'genes_' + version + "/" + ensembl_h + "_" + str(i)
        # Create the directory and its parents if they don't exist
        os.makedirs(new_dir, exist_ok=True)
        ortholgs_output = new_dir + "/orthologs_points.txt"
        statistics_output_h = new_dir + "/human_statistics.txt"
        statistics_output_m = new_dir + "/mouse_statistics.txt"
        junctions = pd.DataFrame(columns=columns_junctions)
        data_df_h = pd.DataFrame(columns=columns_statistics).set_index('ensembl') #statistic about human exons
        data_df_m = pd.DataFrame(columns=columns_statistics).set_index('ensembl')
        #seq of exons
        (exons_seq_h, exon_start_end_h) = get_seq_exons(tra_exons_h, 'cdna', strand_h) #i can't use cds as type for exons seq
        (exons_seq_m, exon_start_end_m) = get_seq_exons(tra_exons_m, 'cdna', strand_m) 
        #align the exons
        df_id = alignment(exons_seq_h, exons_seq_m, new_dir) 
        #find pairs of human-mouse exons and orthologs junctions
        human_id = df_id.iloc[0:len(exons_seq_h),len(exons_seq_h):]
        # Create a new column 'max_column' to store the column names with max values
        human_id['max_column'] = human_id.apply(lambda row: row.index[row >= min_id].tolist(), axis=1)
        human_id['max_value'] = human_id.iloc[:,:-1].max(axis=1)
        mouse_id = df_id.iloc[len(exons_seq_h):, 0:len(exons_seq_h)]
        mouse_id['max_column'] = mouse_id.apply(lambda row: row.index[row >= min_id].tolist(), axis=1)
        mouse_id['max_value'] = mouse_id.iloc[:,:-1].max(axis=1)
        mouse_id_filter = mouse_id[mouse_id['max_value'] >= min_id]
        for exon_m in mouse_id_filter.index:
            for exon_h_from_m_df in mouse_id_filter.loc[exon_m, 'max_column']:
                for exon_m_from_h_df in human_id.loc[exon_h_from_m_df, 'max_column']:
                    human_id_mouse = human_id.loc[exon_h_from_m_df, 'max_value']
                    if exon_m_from_h_df == exon_m and human_id_mouse >= min_id:
                        #get the rank of h and m exons
                        #junctions are the start and end of each exon
                        exons_rank_h = get_rank(canonical_exons_h, exon_h_from_m_df)
                        exons_rank_m = get_rank(canonical_exons_m, exon_m_from_h_df)
                        exon_count_h[exons_rank_h] = 1
                        exon_count_m[exons_rank_m] = 1
                        position_h = exon_start_end_h[exon_h_from_m_df].split("_")
                        position_m = exon_start_end_m[exon_m_from_h_df].split("_")
                        (start_h, end_h) = junction_position(position_h, strand_h, chr_h)
                        (start_m, end_m) = junction_position(position_m, strand_m, chr_m)
                        #'ensembl', 'symbol', 'exon_h', 'rank_h', 'position_h', 'exon_m', 'rank_m', 'position_m'
                        id_human = df_id.loc[exon_m_from_h_df][exon_h_from_m_df]
                        id_mouse = df_id.loc[exon_h_from_m_df][exon_m_from_h_df]
                        average_id = (id_human+id_mouse)/2
                        #the start and end of each orthologous exons are an orthologs point
                        new_row = pd.DataFrame([[ensembl_h,gene_h,exon_h_from_m_df,exons_rank_h,start_h,exon_m_from_h_df,exons_rank_m,start_m,id_human,id_mouse,average_id]], columns=junctions.columns)
                        junctions = pd.concat([junctions, new_row], ignore_index=True)
                        new_row = pd.DataFrame([[ensembl_h,gene_h,exon_h_from_m_df,exons_rank_h,end_h,exon_m_from_h_df,exons_rank_m,end_m,id_human,id_mouse,average_id]], columns=junctions.columns)
                        junctions = pd.concat([junctions, new_row], ignore_index=True)
        if len(junctions)>2: #more than 1 ortholog exons
            junctions.to_csv(ortholgs_output, sep="\t") 
            good_orthologs_h = (human_id['max_value'] >= min_id).sum()
            good_orthologs_m = (mouse_id['max_value'] >= min_id).sum()
            #columns_statistics = ['ensembl','gene','transcripts','exons','exons_with_ortholog','E1','E2']
            data_df_h.loc[ensembl_h] = [gene_h,len(tra_exons_h),len(all_exons_h),good_orthologs_h,exon_count_h['E1'],exon_count_h['E2']]
            data_df_m.loc[ensembl_m] = [gene_m,len(tra_exons_m),len(all_exons_m),good_orthologs_m,exon_count_m['E1'],exon_count_m['E2']]
            data_df_h.to_csv(statistics_output_h, sep="\t")
            data_df_m.to_csv(statistics_output_m, sep="\t")
        else: #less than 2 ortholg exons
            shutil.rmtree(new_dir) #remove the folder
            print(i, gene_h, ensembl_h, ensembl_m, 'less than 2 orthologs exons')
            output = working_dir + "genesWithOneOrthologExon_" + version + ".txt"
            with open(output, "a") as file:
                toWrite = " ".join([str(i), gene_h, ensembl_h, ensembl_m, "\n"])
                file.write(toWrite)                                     
    else: #less than 2 exons
        print(i, gene_h, ensembl_h, ensembl_m, 'less than 2 exons')
        output = working_dir + "genesWithOneOrthologExon_" + version + ".txt"
        with open(output, "a") as file:
            toWrite = " ".join([str(i), gene_h, ensembl_h, ensembl_m, "\n"])
            file.write(toWrite)  
    return
          
if __name__ == "__main__":
      main() 
