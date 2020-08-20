
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

struct Node{
	int label;  // corresponding to an id number of present sample
	double len_above; // time to the parent node; or -1 if it is the root node
	int parent; // index to parent node in the nodeArr array
	int left; // index to left child node in the nodeArr array
	int right;// index to right child node in the nodeArr array
	double time; // time to present
};

int * left_array=NULL;   // temporary array 
int * right_array=NULL;  // temporary array
int countLeft=0, countRight=0; //counter of used items in the array
char * buff=NULL;  // buff for reading a line from msmso file
size_t max_size = 400000000; // max_size of the buffer
struct Node * nodeArr=NULL; // nodeArr for storing parsed data for each tree; reset before each tree parsing
int numNode=0; // number of nodes needed; it is the sam for all trees in a single msms file, thus only needed intialize once.
int nsam=0; // calculated number of present samples; nsam = (numNode + 1)/2
int hom_het = 0;
char * endPtr=NULL; // tempory holder variable for calling functions such strotol 
long true_IBD_sampling_window=0;

/* matrices */
//
// size = nsam * nsam; but only lower trangle were used for efficiency
// each point represent current ibd info for the sample pair (i, j), i is row no., j is col no.
long *ibdStartMatrix=NULL;   // ibd start coord
long *ibdEndMatrix=NULL; // ibd end coord; it is actuall the next base of the ibd's last base, ie end coord not included in the segment
double *ibdTmrcaMatrix=NULL; //  the constant Tmrca for this ibd

/* info of the segment related to the tree being parsed and analyzed. */
long segStart=0; 
long segEnd=0;
double segTmrca=0;
long segLength=0;
long runningLength=0;

/* some constants */
double delta=0;  // tolerance when comparing ibd tmrca with segment tmrca. initialized in the main function
long bp_per_cm = 15000; // related to recomibnation rates
long time_scale_factor=20000;
int min_cm = 2; // cutoff, only IBD longer than min_cm will be printed out.
long min_bp =0; // cutoff, min_bp is calcuated in main function, min_bp = min_cm * bp_per_cm

/* the chrom number of the msmso file */
int chrom=1;  // used second command line parameter to update the chromsome info

/* Files */
FILE * f_map=NULL;
FILE * f_ibd=NULL;
FILE * f_log=NULL;
FILE * f_vcf=NULL;


/* 
 * Recursive function to parse a node info.
 *
 * Key algorithm: 
 * 1. direction: outside to inside, corresponding root to leaves
 * 2. finding the taxon and len_above separator -- ":" or ";", from right to left
 * 3. finding the the seperatro of left child node and right child node -- ',',
 * The critcal condition is that # of right parentheses is 1 more than # of left parentheses
 * between the comma separator and colon separator.
 */
int parse_a_node(char* start, char* end, int * pCount, int parent)
{
	/*
	char *p;
	for (p=start; p<=end; p++) printf("%c", *p);
	printf("\n");
	*/

	// NOTE: every function should have local cur, otherwise the value will
	// be assigned to wrong node
	struct Node * cur=NULL; 

	(*pCount)++; // number of node used
	int index = *pCount -1 ; 
	cur = nodeArr + index;
	cur->parent = parent;
	// find ; or :
	char *pSep=end;
	while(*pSep!=':' && *pSep!=';' && pSep>=start)pSep--;
	if (*pSep == ';') {
		//root sample
		cur -> len_above = 0;
	}
	else if (*pSep == ':')
		//sscanf(pSep+1, "%lf", &(cur->len_above));
		cur->len_above = strtod(pSep+1, &endPtr);
	// find , in outmost ()
	int diff_RLparentheses=0;
	char *pComma=pSep;
	while(!(*pComma==',' && diff_RLparentheses==1)){
		if(*pComma==')') diff_RLparentheses++;
		else if(*pComma=='(') diff_RLparentheses--;
		pComma--;
		if(pComma<=start) break;
	}
	// no , current sample
	if (*pComma != ',') {
		cur->label = atoi(start);
		cur->left = -1;
		cur->right = -1;
	}
	// with , ancestral sample
	else{
		cur->left = parse_a_node(start+1, pComma-1, pCount, index);
		//printf("cur->left; %d\n", cur->left);
		cur->right = parse_a_node(pComma+1, pSep-2, pCount, index);
		//printf("cur->rght; %d\n", cur->right);
		cur->label = -1;
	}
	return index;
}

/*
 * Wrapper of the parse_a_node function to parse the whole tree
 * In addition,
 * 	calculte the time of all ancestral samples
 *
 * 	Key algorithm: 
 * 		use a loop to track down the left child node of each left child node until reaching the present sample
 * 		the total edge length along the path is the time of the ancestral node
 * */
void parse_a_tree(char * tree, long tree_str_len)
{
	//intialize
	int count=0;

	// for(int i=0; i<numNode; i++){
	// 	
	// }
	memset(nodeArr,0, sizeof(struct Node) * numNode);
	// call recursive function to parse the tree
	parse_a_node(tree, tree+tree_str_len-1, &count, -1);

	/* find time for all ancetral node*/
	double t;	
	int i, j;
	for (i=0; i< numNode; i++){
		if(nodeArr[i].left==-1)continue;
		t=0.0;
		for(j=nodeArr[i].left; j!=-1; j=nodeArr[j].left) t += nodeArr[j].len_above;
		nodeArr[i].time=t;
	}
	/* testing 
	for(int i=0; i< count; i++) 
		printf("label: %d, length: %lf, parent: %d, left: %d, right: %d, time: %0.5f\n",
				nodeArr[i].label,
				nodeArr[i].len_above,
				nodeArr[i].parent,
				nodeArr[i].left,
				nodeArr[i].right,
				nodeArr[i].time);
	*/
}

/*
 * A recursive function to find all leaf nodes of the right arm of a ancestral node ancNode (set LEFT) as well as
 * all leaf nodes of the left arm of the ancestral node (set RIGHT).
 *
 * Key algorithm: 
 * 	for any node i in (set LEFT) and any node j in (set RIGHT), 
 * 	tmrca(i, j) is the time of ancNode.
 *
 * */
void find_present_sample(int index, int * array, int * pCount1)
{
	if (nodeArr[index].left==-1){
		array[*pCount1] = nodeArr[index].label;
		(*pCount1)++;
	}
	else {
		find_present_sample(nodeArr[index].left,array, pCount1);
		find_present_sample(nodeArr[index].right, array, pCount1);
	}
}

/*
 * Wrapper of find_present_sample
 * In addition, 
 * 	using segment info and sample info to update the IBD matrices.
 *
 * 	Key algorithm:
 * 	1) for each pair of present sample, from find_present_sample, obtain id1, id2, and trmca
 * 	2) for each tree, we get the segStart, segEnd, info.
 * 	3) if ibdTmrcaMatrix[id1, id2] == trmca, elongate this IBD by setting ibdEndMatrix[id1, id2] = segEnd
 * 	4) otherwise, IBD is not continous with the segment being analysed.
 * 		a) if the IBD is longer than > 2cM right the IBD to stdout
 * 		b) set IBD[id1, id2] to match current segment info
 * */
void calculate_tmrca_and_output_ibd()
{
	double *pDbl=NULL;
	int id1, id2, temp;

	// for each segment or tree, update segment info: length, start and end
	// for each ancester node, update tmrca
	long ibdLength;
	long s, e;

	// for each ancester node i, Find all current sample connected to the
	// left side (left set) this node and find all current samples
	// connected to the right sides (left set). Then any one from left set
	// and another one from right set as a pair will have MRCA == i
	for(int i=0; i<numNode; i++)
	{
		if(nodeArr[i].left == -1) continue;
		countLeft=0, countRight=0;
		for(int i=0; i<numNode; i++) {
			left_array[i]=0;
			right_array[i]=0;
		}
		find_present_sample(nodeArr[i].left, left_array, &countLeft);
		find_present_sample(nodeArr[i].right, right_array, &countRight);

		// update Tmrca for a pair of samples colescet at node i in this segment
		segTmrca = nodeArr[i].time;

		for(int left=0; left < countLeft; left++)
			for (int right=0; right < countRight; right++)
			{
				id1 = left_array[left];
				id2 = right_array[right];
				//only use lower matrix
				if(id1<id2){temp=id1; id1=id2; id2=temp;}
				// printf("id1: %d, id2: %d\n", id1, id2);
				if(fabs(segTmrca-ibdTmrcaMatrix[id1*nsam+id2])<delta) {
					ibdEndMatrix[id1*nsam+id2] = segEnd;
				}

				else{
					e = ibdEndMatrix[id1*nsam+id2];
					s = ibdStartMatrix[id1*nsam+id2];
					ibdLength = e - s;
					if (ibdLength >= min_bp)
						fprintf(f_ibd, "%d\t1\t%d\t1\t%d\t%ld\t%ld\t%0.4f\t3\t%lf\n",
								id1, id2, chrom, s, e, 
								1.0 * ibdLength/bp_per_cm,
								ibdTmrcaMatrix[id1*nsam+id2]
								);
					// update new IBD segment for this paris in the matrics
					ibdStartMatrix[id1*nsam+id2] = segStart;
					ibdEndMatrix[id1*nsam+id2] = segEnd;
					ibdTmrcaMatrix[id1*nsam+id2] = segTmrca;
				}
			}
	}
}

/*
 * Main function:
 *
 * 	1. parse command line parameter, get filename and chrom
 * 	2. loop through the lines contains tree info and segment length info
 * 	3. in the first loop, parse the number of nodes needed by checking
 * 	number of ";" and  ":" separators. Based on numNode, calculate nsam.
 * 	With numNode and nsam node, allocate memory for all the arrays and matrices
 * 	3. in each loop, call the tree parsing function and the function to calculate 
 * 	mrca and output IBD.
 * */
int main(int argc, char * argv[]){

	if(argc<6){
		printf(" Usage. 	macs 50 10000 -t 0.0004 -r 0.0004 -h 1000 -T 2>/dev/null | \\\n"
				"		macs_tree_parser <chromN> <bp_per_cm> <time_scale_factor> <true_IBD_sampling_window> [<hom_het:0>]\n"
				"\n"
				"	chromN:                                            from 1 to 14 for pf; from 1 to 30 for human\n"
				"	bp_per_cm:                                         1000000 for human; 15000 for p.f.\n"
				"	time_scale_factor:                                 4N0 (macs)\n"
				"	true_IBD_sampling_window:                          1 check all trees, 1000, will check tree covers 0, 1000, 2000 \n"
				"		                                           suggestions: human->10000, pf->150; \n" 
				"	hom_het:                                           0 for heterozygous vcf; 1 for homozygous \n\n"); 
		return -1;
	}
	chrom = atoi(argv[1]);
	bp_per_cm=strtol(argv[2], &endPtr,10);
	time_scale_factor = (long)strtod(argv[3], NULL);
	true_IBD_sampling_window = strtol(argv[4], NULL, 10);
	if(true_IBD_sampling_window==0) {fprintf(stderr, "true_IBD_sampling_window parameter error!\n"); exit(1);}
	if(argc>=6)
          hom_het = atoi(argv[5]);

        delta=0.5/time_scale_factor;  // tolerance when comparing ibd tmrca with segment tmrca.
	char filename[10];
	sprintf(filename, "%d.ibd", chrom); f_ibd = fopen(filename, "w");
	sprintf(filename, "%d.log", chrom); f_log = fopen(filename, "w");
	sprintf(filename, "%d.map", chrom); f_map = fopen(filename, "w");
	sprintf(filename, "%d.vcf", chrom); f_vcf = fopen(filename, "w");

	size_t lineSize=max_size;
	char * tree = NULL;
	long tree_str_len = 0;
	char * token= NULL;
	char * pChar=NULL;
	long position=1, pre_position=1, last_print_pos=1;
	long segno=0;
	long nsites=0;

	segEnd=0, segStart=0;
	min_bp = min_cm * bp_per_cm;
	buff = (char*) malloc(sizeof(char)*max_size);

	long counter=0;
	while (getline(&buff, &lineSize, stdin) >= 0) {
		// check if line is cmdline and determine length info and sample size
		if (strncmp(buff, "COMMAND:", 8) == 0) {
			for(pChar=buff; *pChar!='m'; ++pChar){};
			token = strtok(pChar, " "); // macs
			token = strtok(NULL, " "); // nsam
			nsam = (int)strtod(token, NULL);
			numNode = 2 * nsam - 1;
			token = strtok(NULL, " "); // nsites
			nsites=(long)strtod(token, NULL);

			// alloc memory for nodeArr, initialize numNode
			if (nodeArr == NULL) {
				nodeArr = (struct Node *)malloc(sizeof(struct Node) * numNode);
				left_array = (int *)malloc(sizeof(int) * numNode);
				right_array = (int *)malloc(sizeof(int) * numNode);
				// matrices
				ibdStartMatrix = (long *)calloc(nsam * nsam, sizeof(long));
				ibdEndMatrix = (long *)calloc(nsam * nsam, sizeof(long));
				ibdTmrcaMatrix = (double *)calloc(nsam * nsam, sizeof(double));

				for (int i = 0; i < nsam; i++)
					for (int j = 0; j < i; j++) {
						// using calloc or memset can not correctly intialize arrays to
						// 0 need to explict initialize it using loop
						ibdStartMatrix[i * nsam + j] = 0;
						ibdEndMatrix[i * nsam + j] = 0;
						ibdTmrcaMatrix[i * nsam + j] = 0;
					}
			}
		}
		// check if segsites line
		else if (strncmp(buff, "SITE:", 5) == 0) {
			// find nseg
			token = strtok(buff, " \t"); // SITE: 
			token = strtok(NULL, " \t"); // segno
			token = strtok(NULL, " \t"); // poistion
			pre_position = position;
			position = (long)(strtod(token, NULL) * nsites);
			strtok(NULL, " \t"); // the 3rd num; do know what it is.
			pChar = strtok(NULL, " \t"); // start of 01
			if (position <= pre_position) position = pre_position + 1;

			// map file ( will not print every segsite, only print one position every 0.01cM
			if(position > last_print_pos + bp_per_cm / 100)
			{
				fprintf(f_map, "%d\tsnp_%d_%ld\t%lf\t%ld\n", 
						chrom, 
						chrom, position, 
						position * 1.0 /bp_per_cm,
						position );
				last_print_pos = position;
			}

			
			// ------ write vcf file
			// header 
			// 	 ## lines + #CHROM line
			if(segno == 0)
			{
				fprintf(f_vcf,
						"##fileformat=VCFv4.2\n" 
						"##source=macs_tree_parser by B.G.\n"
						"##contig=<ID=%d,length=%ld>\n"
						"##INFO=<ID=PR,Number=1,Type=Flag,Description=\"coverted from ms format data\">\n"
						"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
						chrom, nsites);
				//       ## sample names
				if (hom_het == 0) // het
					for (long i = 0; i < nsam / 2; i++) fprintf(f_vcf, "\t%ld", i);
				else // hom
					for (long i = 0; i < nsam; i++) fprintf(f_vcf, "\t%ld", i);

				fprintf(f_vcf, "\n");
			}

			// Body
			// 	find two different alleles
			char ATGC[] = "ATGC";
			int ref, alt;
			ref = rand() % 4;
			alt = ref;
			while(alt == ref) alt = rand() % 4;
			//	printf intial columns
			fprintf(f_vcf, "%d\t%ld\tsnp_%d_%ld\t%c\t%c\t.\t.\tPR\tGT",
					chrom,
					position,
					chrom, position,
					ATGC[ref], ATGC[alt]);
			//	printf allele columns
			if(hom_het==0) // het
				for(int col=0; col< nsam; col+=2) fprintf(f_vcf, "\t%c|%c", 
					pChar[col], pChar[col+1]);
			else // homo
				for(int col=0; col< nsam; ++col) fprintf(f_vcf, "\t%c|%c", 
					pChar[col], pChar[col]);
			fprintf(f_vcf,"\n");
			++segno;
		}
		// --------------
		// for lines with tree info
		else if (strncmp(buff, "NEWICK_TREE:", 12) == 0)
		{
			// Get seg length
			for(pChar=buff; *pChar!='['; ++pChar){}; ++pChar; // move to seg length position
			segLength = strtol(pChar, NULL, 10); 
			segStart= segEnd;
			segEnd += segLength;

			// if a mulplte of true_IBD_sampling_window is not on the segment, skip parsing the tree
			// using browning method
			if(segEnd/true_IBD_sampling_window == segStart/true_IBD_sampling_window) {
				continue;
			}

			for (tree = pChar; *tree != '('; tree++){} ; // find the first (
			tree_str_len = strlen(tree);

			// parse the tree using recursive function
			parse_a_tree(tree, tree_str_len);

			// find ancetral nodes 's sample pairs whose mrca is this node
			// then update ibd matrix and output ibd
			calculate_tmrca_and_output_ibd();

		}
		// logging 
		counter++;
		fprintf(f_log, "nsam: %d, TreeCounter: %ld, segStart: %ld\n", nsam,
				counter, segStart);

		// reset line max size for next lines
		lineSize = max_size;
	}
	fclose(f_ibd);
	fclose(f_map);
	fclose(f_log);

	if (nodeArr != NULL) {
		free(nodeArr); nodeArr=NULL;
		free(left_array); left_array=NULL;
		free(right_array); right_array=NULL;
		free(ibdStartMatrix); ibdStartMatrix=NULL;
		free(ibdEndMatrix); ibdEndMatrix=NULL;
		free(ibdTmrcaMatrix); ibdTmrcaMatrix=NULL;
	}
	pChar=NULL; endPtr=NULL;

	if(buff!=NULL){
		free(buff);
		buff=NULL;
	}

	return 0;
}
