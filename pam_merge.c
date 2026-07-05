#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>
#include <sys/queue.h>
#include <dotgeno.h>

#define PAM_STR_BUF_EXTRA 6

// important global variables
bool HASH_CHECK = true;
bool IS_VERBOSE = false;

typedef struct {
	snp_data* elems;
	size_t length;
} snp_objs;

typedef struct {
	ind_data* elems;
	size_t length;
} ind_objs;

typedef struct {
	pam_file_reader* elems;
	size_t length;
} pam_objs;

typedef struct {
	struct idx_head* ref_idx;
	struct idx_head* elm_idx;  // come up with better name
} idx_pair;

typedef struct {
	idx_pair** ip;
	size_t length;
	bool is_intersected;
} idx_intersect;

typedef enum {
	MATCH,
	FLIP,
	MISMATCH
} snp_cmp_res;

typedef struct {
	char** snp;
	char** ind;
	char** pam;
	size_t length;
} combined_files;

snp_objs init_snp_objs(size_t length) {
	snp_objs snps;
	snps.length = length;
	snps.elems = (snp_data*)malloc(sizeof(snp_data) * snps.length);
	return snps;
}

ind_objs init_ind_objs(size_t length) {
	ind_objs inds;
	inds.length = length;
	inds.elems = (ind_data*)malloc(sizeof(ind_data) * inds.length);
	return inds;
}

pam_objs init_pam_objs(size_t length) {
	pam_objs pams;
	pams.length = length;
	pams.elems = (pam_file_reader*)malloc(sizeof(pam_file_reader) * pams.length);
	return pams;
}

idx_intersect init_idx_intersect(snp_objs* snps) {
	idx_intersect isct;
	isct.length = snps->length - 1;
	isct.is_intersected = false;
	isct.ip = (idx_pair**)malloc(sizeof(idx_pair*) * isct.length);
	return isct;
}

idx_pair* get_idx_pair(snp_data* ref, snp_data* cmp) {
	idx_pair* pair =(idx_pair*) malloc(sizeof(idx_pair));
	pair->ref_idx = (struct idx_head*)malloc(sizeof(struct idx_head));
	pair->elm_idx = (struct idx_head*)malloc(sizeof(struct idx_head));

	STAILQ_INIT(pair->ref_idx);
	STAILQ_INIT(pair->elm_idx);
	intersect_snp_data(ref, cmp, pair->ref_idx, pair->elm_idx);
	return pair;
}

size_t get_max_index(struct idx_node** arr, size_t length) {
	size_t max = 0;
	size_t max_idx = 0;
	for(size_t i = 0; i < length; i++) {
		if(arr[i]->idx >= max) {
			max = arr[i]->idx;
			max_idx = i;
		}
	}
	return max_idx;
}

bool all_equal(struct idx_node** arr, size_t length) {
	size_t first = arr[0]->idx;
	for(size_t i = 1; i < length; i++) {
		if(first != arr[i]->idx) {
			return false;
		}
	}
	return true; 
}

void print_isct(idx_intersect* isct) {
	for(size_t i = 0; i < isct->length; i++) {
		printf("Index %u:\n", i);
		printf("\tIDX:");
		struct idx_node* in;
		STAILQ_FOREACH(in, isct->ip[i]->ref_idx, nodes) {
			printf(" %u", in->idx);
		}
		printf("\n");
	}
}

void make_intersection(idx_intersect* isct) {
	// store the current nodes in array
	struct idx_node** cur_nodes_ref = (struct idx_node**)malloc(isct->length * sizeof(struct idx_node*));
	struct idx_node** cur_nodes_elm = (struct idx_node**)malloc(isct->length * sizeof(struct idx_node*));
	for(size_t i = 0; i < isct->length; i++) {
		cur_nodes_ref[i] = STAILQ_FIRST(isct->ip[i]->ref_idx);
		cur_nodes_elm[i] = STAILQ_FIRST(isct->ip[i]->elm_idx);
	}

	print_isct(isct);
	while(1) {
		bool end_while = false;  // change value to exit while loop
		// if all the nodes are equal then move on to the next elements
		if(all_equal(cur_nodes_ref, isct->length)) {
			for(size_t i = 0; i < isct->length; i++) {
				struct idx_node* next_node_ref = STAILQ_NEXT(cur_nodes_ref[i], nodes);
				struct idx_node* next_node_elm = STAILQ_NEXT(cur_nodes_elm[i], nodes);
				// if this is the last element then remove all elements after this current one from the lists
				if(next_node_ref == NULL) {
					for(size_t j = i + 1; j < isct->length; j++) {
						struct idx_node* nxt_r = STAILQ_NEXT(cur_nodes_ref[j], nodes);
						struct idx_node* nxt_e = STAILQ_NEXT(cur_nodes_elm[j], nodes);
						while(nxt_r) {
							STAILQ_REMOVE(isct->ip[j]->ref_idx, nxt_r, idx_node, nodes);
							STAILQ_REMOVE(isct->ip[j]->elm_idx, nxt_e, idx_node, nodes);
							free(nxt_r);
							free(nxt_e);
							nxt_r = STAILQ_NEXT(cur_nodes_ref[j], nodes);
							nxt_e = STAILQ_NEXT(cur_nodes_elm[j], nodes);
						}
					}
					for(size_t j = 0; j < i; j++) {
						struct idx_node* nxt_r = cur_nodes_ref[j];
						struct idx_node* nxt_e = cur_nodes_elm[j];
						while(nxt_r) {
							STAILQ_REMOVE(isct->ip[j]->ref_idx, nxt_r, idx_node, nodes);
							STAILQ_REMOVE(isct->ip[j]->elm_idx, nxt_e, idx_node, nodes);
							cur_nodes_ref[j] = nxt_r;
							cur_nodes_elm[j] = nxt_e;
							nxt_r = STAILQ_NEXT(nxt_r, nodes);
							nxt_e = STAILQ_NEXT(nxt_e, nodes);
							free(cur_nodes_ref[j]);
							free(cur_nodes_elm[j]);
						}
					}
					end_while = true;
					break;
				}
				cur_nodes_ref[i] = next_node_ref;
				cur_nodes_elm[i] = next_node_elm;
			}
		// if not all the nodes are equal then remove all current elements EXCEPT the 
		// maximum ones and move on to the next elements
		} else {
			size_t max_i = get_max_index(cur_nodes_ref, isct->length);
			for(size_t i = 0; i < isct->length; i++) {
				// ignore max values
				if(i == max_i) { continue; }
				if(cur_nodes_ref[i]->idx == cur_nodes_ref[max_i]->idx) { continue; }
				// remove and move to the next one
				struct idx_node* next_node_ref = STAILQ_NEXT(cur_nodes_ref[i], nodes);
				struct idx_node* next_node_elm = STAILQ_NEXT(cur_nodes_elm[i], nodes);
				STAILQ_REMOVE(isct->ip[i]->ref_idx, cur_nodes_ref[i], idx_node, nodes);
				STAILQ_REMOVE(isct->ip[i]->elm_idx, cur_nodes_elm[i], idx_node, nodes);
				free(cur_nodes_ref[i]);
				free(cur_nodes_elm[i]);
				// if one of the next nodes is empty list and free elements for all of them, there is no intersection
				if(next_node_ref == NULL) {
					for(int j = 0; j < isct->length; j++) {
						free_idx_list(isct->ip[j]->ref_idx);
						free_idx_list(isct->ip[j]->elm_idx);
					}
					end_while = true;
					break;
				}
				cur_nodes_ref[i] = next_node_ref;
				cur_nodes_elm[i] = next_node_elm;
			}
		}
		print_isct(isct);
		if(end_while) {
			free(cur_nodes_ref);
			free(cur_nodes_elm);
			break;
		}
	}
	isct->is_intersected = true;
}

ind_data append_ind_objs(ind_objs* inds) {
	ind_data ind_total = inds->elems[0];
	for(size_t i = 1; i < inds->length; i++) {
		ind_data out_buf;
		append_ind_data(&ind_total, &inds->elems[i], &out_buf);
		if(i > 1) { free_ind_data(&ind_total); }
		ind_total = out_buf;
	}
	return ind_total;
}

void free_snp_objs(snp_objs* snps) {
	for(size_t i = 0; i < snps->length; i++) {
		free_snp_data(&snps->elems[i]);
	}
	free(snps->elems);
}

void free_ind_objs(ind_objs* inds) {
	for(size_t i = 0; i < inds->length; i++) {
		free_ind_data(&inds->elems[i]);
	}
	free(inds->elems);
}

void free_pam_objs(pam_objs* pams) {
	for(size_t i = 0; i < pams->length; i++) {
		close_pam_file_reader(&pams->elems[i]);
	}
	free(pams->elems);
}

void free_idx_pair(idx_pair* pair) {
	free_idx_list(pair->ref_idx);
	free_idx_list(pair->elm_idx);
	free(pair->ref_idx);
	free(pair->elm_idx);
	free(pair);
}

void free_idx_intersect(idx_intersect* isct) {
	for(size_t i = 0; i < isct->length; i++) {
		free_idx_pair(isct->ip[i]);
	}
	free(isct->ip);
}

short get_snp_data_from_idx_intersect(snp_data* snp_in, snp_data* snp_out, idx_intersect* isct) {
	if(!isct->is_intersected) {
		fprintf(stderr, "ERROR: SNP files have not been intersected.");
		exit(EXIT_FAILURE);
	}
	return filter_snp_data(snp_in, snp_out, isct->ip[0]->ref_idx);
}


// switch to a generic method return snp_objs
snp_objs filter_snp_objs(snp_objs* snps) {
	idx_intersect isct = init_idx_intersect(snps);
	for(size_t i = 0; i < isct.length; i++) {
		isct.ip[i] = get_idx_pair(&snps->elems[0], &snps->elems[i+1]);
	}
	make_intersection(&isct);
	snp_objs s_out;
	s_out = init_snp_objs(snps->length);
	filter_snp_data(&snps->elems[0], &s_out.elems[0], isct.ip[0]->ref_idx);
	for(size_t i = 1; i < s_out.length; i++) {
		filter_snp_data(&snps->elems[i], &s_out.elems[i], isct.ip[i-1]->elm_idx);
	}
	free_idx_intersect(&isct);
	return s_out;
}

snp_cmp_res snp_cmp(snp_data* snp1, snp_data* snp2, size_t idx1, size_t idx2) {
	if((strcmp(snp1->ref[idx1], snp2->ref[idx2]) == 0) && (strcmp(snp1->alt[idx1], snp2->alt[idx2]) == 0)) {
		return MATCH; 
	} else if((strcmp(snp1->ref[idx1], snp2->alt[idx2]) == 0) && (strcmp(snp1->alt[idx1], snp2->ref[idx2]) == 0)) {
		return FLIP;
	}
	else {
		return MISMATCH;
	}
}

void write_combined_pam(pam_objs* pams, snp_objs* snps, ind_objs* inds, char* out_prefix) {
	snp_data snp_ref = snps->elems[0];
	ind_data ind_total = append_ind_objs(inds);
	char* buf = (char*)malloc(sizeof(char) * (strlen(out_prefix) + PAM_STR_BUF_EXTRA));
	sprintf(buf, "%s.snp", out_prefix);
	write_snp_data(&snp_ref, buf);
	sprintf(buf, "%s.ind", out_prefix);
	write_ind_data(&ind_total, buf);
	sprintf(buf, "%s.geno", out_prefix);	
	pam_file_writer paw = pam_file_writer_init(buf, &snp_ref, &ind_total);
	free(buf);
	if(!(pams->length == snps->length) && (snps->length == inds->length)) {
		fprintf(stderr, "ERROR: inconsistent lengths for pams, snps, and inds.\n");
		exit(EXIT_FAILURE);
	}
	if(!paw.is_open) {
		fprintf(stderr, "ERROR: cannot merge closed pam_file_writer.\n");
		exit(EXIT_FAILURE);
	}
	if(paw.hdr_written) {
		fprintf(stderr, "ERROR: header of PACKEDANCESTRYMAP already written.\n");
		exit(EXIT_FAILURE);
	}
	// read pam headers and check hashes if specified
	for(size_t i = 0; i < pams->length; i++) {
		hdr_data hdr = read_pam_header(&pams->elems[i]);
		if(!HASH_CHECK) { continue; }
		if(hdr.ind_hash != inds->elems[i].hash) {
			fprintf(stderr, "Incorrect hash in one of the ind files. ind file may have been modiied after PACKEDANCESTRYMAP file creation.\n");
			exit(EXIT_FAILURE);
		}
		if(hdr.snp_hash != snps->elems[i].hash) {
			fprintf(stderr, "Incorrect hash in one of the snp files. snp file may have been modiied after PACKEDANCESTRYMAP file creation.\n");
			exit(EXIT_FAILURE);
		}
	}

	write_pam_header(&paw, &snp_ref, &ind_total);
	for(size_t i = 0; i < snp_ref.length; i++) {
		uint8_t* record = (uint8_t*)malloc(sizeof(uint8_t) * ind_total.length);
		size_t adder = 0;
		for(size_t j = 0; j < pams->length; j++) {
			short ret = goto_var_pam(&pams->elems[j], &snps->elems[j], snp_ref.var_id[i]);
			if(ret == 0) {
				snp_cmp_res cmp_res = snp_cmp(&snp_ref, &snps->elems[j], i, pams->elems[j].idx);
				uint8_t* rcd = read_pam_record(&pams->elems[j]);
				switch (cmp_res) {
					case MATCH:
						for(size_t k = 0; k < inds->elems[j].length; k++) {
							record[k + adder] = rcd[k];
						}
						break;
					case FLIP:
						for(size_t k = 0; k < inds->elems[j].length; k++) {
							record[k + adder] = ((rcd[k] != NAN_VAL)*(2-rcd[k])) + ((rcd[k] == NAN_VAL)*NAN_VAL);
						}
						break;
					case MISMATCH:
						for(size_t k = 0; k < inds->elems[j].length; k++) {
							record[k + adder] = NAN_VAL;
						}
				}
				free(rcd);
			} else if(ret == -1) {
				for(size_t k = 0; k < inds->elems[j].length; k++) {
					record[k + adder] = NAN_VAL;
				}
			}
			adder += inds->elems[j].length;
		}
		write_pam_record(&paw, record);
		free(record);
	}
	close_pam_file_writer(&paw);
	free_ind_data(&ind_total);
}

char** get_comma_separated_values(char* str, size_t* n_chars) {
	size_t str_len = strlen(str);
	size_t nc_buf = str_len > 0;
	// get num chars
	for(size_t i = 0; i < str_len; i++) {
		if(str[i] == ',') {
			nc_buf += 1;
		}
	}
	size_t* col_lens = (size_t*)malloc(nc_buf * sizeof(size_t)); 
	size_t cur_len = 0;
	size_t cur_col = 0;
	for(size_t i = 0; i < str_len; i++) {
		if(!(str[i] == ',')) {
			cur_len += 1;
		} else {
			col_lens[cur_col] = cur_len;
			cur_len = 0;
			cur_col += 1;
		}
	}
	col_lens[cur_col] = cur_len;
	char** out_arr = (char**)malloc(nc_buf * sizeof(char*));
	for(size_t i = 0; i < nc_buf; i++) {
		out_arr[i] = (char*)malloc((col_lens[i] + 1) * sizeof(char));
		out_arr[i][col_lens[i]] = '\0';
	}
	size_t subber = 0;
	cur_col = 0;
	for(size_t i = 0; i < str_len; i++) {
		if((str[i] == ',')) {
			cur_col += 1;
			subber = i + 1;
		} else {
			out_arr[cur_col][i - subber] = str[i];
		}
	}
	free(col_lens);
	*n_chars = nc_buf;
	return out_arr;
}

void free_gsf(char** arr, size_t len) {
	for(size_t i = 0; i < len; i++) {
		free(arr[i]);
	}
	free(arr);
}

void free_combined_files(combined_files* file_info) {
	free_gsf(file_info->snp, file_info->length);
	free_gsf(file_info->ind, file_info->length);
	free_gsf(file_info->pam, file_info->length);
}

void print_help() {
	printf("Future help message here\n");
}

int main(int argc, char* argv[]) {
	bool intersect_snps = false;
	combined_files file_info;
	file_info.length = 0;
	file_info.snp = NULL;
	file_info.ind = NULL;
	file_info.pam = NULL;
	char* out_prefix = NULL;

	static struct option long_options[] = {
		{"help",           no_argument,       NULL, 'h'},
		{"prefixes",       required_argument, NULL, 'p'},
		{"pams",           required_argument, NULL, 'P'},
		{"snps",           required_argument, NULL, 's'},
		{"inds",           required_argument, NULL, 'i'},
		{"out",            required_argument, NULL, 'o'},
		{"ignore-hash",    no_argument,       0,    300},
		{"intersect-snps", no_argument,       0,    400},
		{"verbose",        no_argument,       0,    500},
		{0,                0,                 0,      0}
	};

	while(1) {
		int c = getopt_long(argc, argv, "hp:P:s:i:o:", long_options, NULL);
		if(c == -1) { break; }
		size_t n;
		char** gsf;
		switch(c) {
			case 'h':
				print_help();
				exit(EXIT_SUCCESS);
			case 'p':
				if(file_info.pam) {
					fprintf(stderr, "ERROR: PACKEDANCESTRYMAP files already provided!");
					exit(EXIT_FAILURE);
				}
				if(file_info.snp) {
					fprintf(stderr, "ERROR: SNP files already provided!");
					exit(EXIT_FAILURE);
				}
				if(file_info.ind) {
					fprintf(stderr, "ERROR: Individual files already provided!");
					exit(EXIT_FAILURE);
				}
				gsf =  get_comma_separated_values(optarg, &n);
				file_info.length = n;
				file_info.snp = (char**)malloc(sizeof(char*) * file_info.length);
				file_info.ind = (char**)malloc(sizeof(char*) * file_info.length);
				file_info.pam = (char**)malloc(sizeof(char*) * file_info.length);
				for(size_t i = 0; i < file_info.length; i++) {
					char* str_snp = (char*)malloc(sizeof(char) * (strlen(gsf[i]) + PAM_STR_BUF_EXTRA));
					char* str_ind = (char*)malloc(sizeof(char) * (strlen(gsf[i]) + PAM_STR_BUF_EXTRA));
					char* str_pam = (char*)malloc(sizeof(char) * (strlen(gsf[i]) + PAM_STR_BUF_EXTRA));
					sprintf(str_snp, "%s.snp", gsf[i]);
					sprintf(str_ind, "%s.ind", gsf[i]);
					sprintf(str_pam, "%s.geno", gsf[i]);
					file_info.snp[i] = str_snp;
					file_info.ind[i] = str_ind;
					file_info.pam[i] = str_pam;
				}
				free_gsf(gsf, n);
				break;
			case 'P':
				if(file_info.pam) {
					fprintf(stderr, "ERROR: Multiple sets of PACKEDANCESTRYMAP files provided!\n");
					exit(EXIT_FAILURE);
				}
				gsf =  get_comma_separated_values(optarg, &n);
				if(file_info.length == 0) {
					file_info.length = n;
				} else {
					if(n != file_info.length) {
						fprintf(stderr, "ERROR: Incompatible number of PACKEDANCESTRYMAP files and SNP and/or individual files!\n");
						exit(EXIT_FAILURE);
					}
				}
				file_info.pam = gsf;
				break;
			case 's':
				if(file_info.ind) {
					fprintf(stderr, "ERROR: Multiple sets of PACKEDANCESTRYMAP files provided!\n");
					exit(EXIT_FAILURE);
				}
				gsf =  get_comma_separated_values(optarg, &n);
				if(file_info.length == 0) {
					file_info.length = n;
				} else {
					if(n != file_info.length) {
						fprintf(stderr, "ERROR: Incompatible number of individual files and PACKEDANCESTRYMAP and/or SNP files!\n");
						exit(EXIT_FAILURE);
					}
				}
				file_info.ind = gsf;
				break;
			case 'i':
				if(file_info.snp) {
					fprintf(stderr, "ERROR: Multiple sets of PACKEDANCESTRYMAP files provided!\n");
					exit(EXIT_FAILURE);
				}
				gsf =  get_comma_separated_values(optarg, &n);
				if(file_info.length == 0) {
					file_info.length = n;
				} else {
					if(n != file_info.length) {
						fprintf(stderr, "ERROR: Incompatible number of SNP files and PACKEDANCESTRYMAP and/or individual files!\n");
						exit(EXIT_FAILURE);
					}
				}
				file_info.snp = gsf;
				break;
			case 'o':
				if(out_prefix) {
					fprintf(stderr, "ERROR: output prefix of '%s' has already been specified!\n", out_prefix);
					exit(EXIT_FAILURE);
				} else {
					out_prefix = optarg;
				}
				break;
			case 300:
				HASH_CHECK = false;
				break;
			case 400:
				intersect_snps = true;
				break;
			case 500:
				IS_VERBOSE = true;
				break;
			case '?':
				break;
			default:
				abort();
		}
	}

	if(file_info.length == 0) {
		fprintf(stderr, "ERROR: no input files provided!\n");
		exit(EXIT_FAILURE);
	}

	if(out_prefix == NULL) {
		fprintf(stderr, "ERROR: no output file prefix provided!\n");
		exit(EXIT_FAILURE);
	}

	// ADD MORE STUFF HERE 
	free_combined_files(&file_info);
	/*
	snp_objs snps = init_snp_objs(argc - 1);
	ind_objs inds = init_ind_objs(argc - 1);
	pam_objs pams = init_pam_objs(argc - 1);
	for(size_t i = 0; i < snps.length; i++) {
		char* buf = (char*)malloc(sizeof(char) * (strlen(argv[i+1]) + PAM_STR_BUF_EXTRA));
		sprintf(buf, "%s.snp", argv[i+1]);
		snps.elems[i] = read_snp_file(buf);
		sprintf(buf, "%s.ind", argv[i+1]);
		inds.elems[i] = read_ind_file(buf);
		sprintf(buf, "%s.geno", argv[i+1]);
		pams.elems[i] = pam_file_reader_init(buf, &snps.elems[i], &inds.elems[i]);
		free(buf);
	}
	write_combined_pam(&pams, &snps, &inds, "out_test");
	free_snp_objs(&snps);
	free_ind_objs(&inds);
	free_pam_objs(&pams);
	*/
}
