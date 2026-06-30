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
	struct idx_head* ref_idx;
	struct idx_head* elm_idx;  // come up with better name
} idx_pair;

typedef struct {
	idx_pair** ip;
	size_t length;
} idx_intersect;

snp_objs init_snp_objs(size_t length) {
	snp_objs snps;
	snps.length = length;
	snps.elems = (snp_data*)malloc(sizeof(snp_data) * snps.length);
	return snps;
}

idx_intersect init_idx_intersect(snp_objs* snps) {
	idx_intersect isct;
	isct.length = snps->length - 1;
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
							nxt_r = STAILQ_NEXT(cur_nodes_ref[j], nodes);
							nxt_e = STAILQ_NEXT(cur_nodes_elm[j], nodes);
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
}

void free_snp_objs(snp_objs* s_objs) {
	for(size_t i = 0; i < s_objs->length; i++) {
		free_snp_data(&s_objs->elems[i]);
	}
	free(s_objs->elems);
}

void free_idx_pair(idx_pair* pair) {
	free_idx_list(pair->ref_idx);
	free_idx_list(pair->elm_idx);
	free(pair);
}

void free_idx_intersect(idx_intersect* isct) {
	for(size_t i = 0; i < isct->length; i++) {
		free_idx_pair(isct->ip[i]);
	}
	free(isct->ip);
}

int main(int argc, char* argv[]) {
	snp_objs snps = init_snp_objs(argc - 1);
	for(size_t i = 0; i < snps.length; i++) {
		snps.elems[i] = read_snp_file(argv[i+1]);
	}
	idx_intersect isct = init_idx_intersect(&snps);
	for(size_t i = 0; i < isct.length; i++) {
		isct.ip[i] = get_idx_pair(&snps.elems[0], &snps.elems[i+1]);
	}
	make_intersection(&isct);
	free_snp_objs(&snps);
	free_idx_intersect(&isct);
}
