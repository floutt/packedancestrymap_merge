#include <stdio.h>
#include <stdlib.h>
#include <sys/queue.h>
#include <stdbool.h>
#include <dotgeno.h>

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

void master_merge(pam_objs* pams, snp_objs* snps, ind_objs* inds, pam_file_writer* paw) {
	snp_data snp_ref = snps->elems[0];
	ind_data ind_total = append_ind_objs(inds);
	if(!(pams->length == snps->length) && (snps->length == inds->length)) {
		fprintf(stderr, "ERROR: inconsistent lengths for pams, snps, and inds.\n");
		exit(EXIT_FAILURE);
	}
	if(!paw->is_open) {
		fprintf(stderr, "ERROR: cannot merge closed pam_file_writer.\n");
		exit(EXIT_FAILURE);
	}
	if(paw->hdr_written) {
		fprintf(stderr, "ERROR: header of PACKEDANCESTRYMAP already written.\n");
		exit(EXIT_FAILURE);
	}
	write_pam_header(paw, &snp_ref, &ind_total);
	for(size_t i = 0; i < snp_ref.length; i++) {
		uint8_t* record = (uint8_t*)malloc(sizeof(uint8_t) * ind_total.length);
		size_t adder = 0;
		for(size_t j = 0; j < pams->length; j++) {
			short ret = goto_var_pam(&pams->elems[j], &snps->elems[j], snp_ref.var_id[i]);
			if(ret == 0) {
				// TODO: check for SNP flips, mismatches, etc.
				uint8_t* rcd = read_pam_record(&pams->elems[j]);
				for(size_t k = 0; k < inds->elems[j].length; k++) {
					record[k + adder] = rcd[k];
					adder += inds->elems[j].length;
				}
				free(rcd);
			} else if(ret == -1) {
				for(size_t k = 0; k < inds->elems[j].length; k++) {
					record[k + adder] = NAN_VAL;
					adder += inds->elems[j].length;
				}
			}
		}
	}
}

int main(int argc, char* argv[]) {
	snp_objs snps = init_snp_objs(argc - 1);
	ind_objs inds = init_ind_objs(argc - 1);
	pam_objs pams = init_pam_objs(argc - 1);
	for(size_t i = 0; i < snps.length; i++) {
		char* buf = (char*)malloc(sizeof(char) * (strlen(argv[i+1]) + 6));
		sprintf(buf, "%s.snp", argv[i+1]);
		snps.elems[i] = read_snp_file(buf);
		sprintf(buf, "%s.ind", argv[i+1]);
		inds.elems[i] = read_ind_file(buf);
		sprintf(buf, "%s.geno", argv[i+1]);
		pams.elems[i] = pam_file_reader_init(buf, &snps.elems[i], &inds.elems[i]);
	}
	
	free_snp_objs(&snps);
	free_ind_objs(&inds);
	free_pam_objs(&pams);
	/*
	snp_objs snps_filtered = filter_snp_objs(&snps);
	free_snp_objs(&snps_filtered);
	free_snp_objs(&snps);
	*/

	/*
	ind_objs inds = init_ind_objs(argc - 1);
	for(int i = 1; i < argc; i++) {
		inds.elems[i-1] = read_ind_file(argv[i]);
	}
	ind_data ind_total = append_ind_objs(&inds);
	write_ind_data(&ind_total, "/dev/stdout");
	free_ind_data(&ind_total);
	free_ind_objs(&inds);
	*/
}
