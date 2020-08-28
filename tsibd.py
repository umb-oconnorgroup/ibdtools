import numpy as np
import msprime
import sys
import pandas as pd


class ts_ibd:
    def __init__(self, tree_sequence, meta_info):

        # make matrics
        num_samples = tree_sequence.num_samples
        self.num_samples = num_samples
        self.tmrca_matrix = np.zeros(shape=(num_samples, num_samples))
        self.breakpos_matrix = np.zeros(shape=(num_samples, num_samples))

        # initialize counter
        self.longest = 0
        self.num_ibd = 0
        self.num_pairs_updated_incremental = 0
        self.num_pairs_updated_complete = 0

        # iterators
        self.trees = tree_sequence.trees()
        self.diffs = tree_sequence.edge_diffs()

        self.tree_sequence = tree_sequence

        # meta info
        self.chrom = meta_info["chrom"]
        self.bp_per_cM = meta_info["bp_per_cM"]
        self.ibd_threshold_in_bp = meta_info["ibd_threshold_in_bp"]

        # running tree and diff, skip the first tree in the for loop
        self.tree = next(self.trees)
        self.diff = next(self.diffs)

        # IBD lists
        self.id1_list = []
        self.id2_list = []
        self.start_list = []
        self.end_list = []
        self.length_list = []
        self.tmrca_list = []

    def initialize_tmrca_matrix_from_first_tree(self):

        first_tree = self.tree_sequence.first()
        self._update_tmrca_matrix_due_to_new_tree(first_tree)

    def iterate_trees_to_update_matrixes(self):

        for i in range(1, tree_sequence.num_trees):
            # get the edge_diff and tree
            self.diff = next(self.diffs)
            self.tree = next(self.trees)

            # if int(interval[0]) / 10000 == int(interval[1]) / 10000:
            #     continue

            print("tree: %d, position: %f, longest: %f, num_ibd: %f" %
                  (i, self.tree.interval[1], self.longest, self.num_ibd),
                  file=sys.stderr)

            M, N, RN, pat = self._find_edge_diff_pattern(self.diff, self.tree)

            self.num_pairs_updated_incremental += \
                self._update_matries_due_to_new_leaves(M, RN)

        self.output_last_ibds_for_all_pairs()

        num_pairs = self.num_samples * (self.num_samples - 1) / 2
        num_pairs_updated_complete = num_pairs * self.tree_sequence.num_trees

        print("", file=sys.stderr)
        print("-----------------", file=sys.stderr)
        print("sample_size: %d" % self.num_samples, file=sys.stderr)
        print("num pairs updated complete   : %d" %
              num_pairs_updated_complete,
              file=sys.stderr)
        print("num pairs updated incremental: %d" %
              self.num_pairs_updated_incremental,
              file=sys.stderr)
        print("incremental to complete ratio: %g" %
              (self.num_pairs_updated_incremental /
               num_pairs_updated_complete),
              file=sys.stderr)

    def _find_tree_leaves_of_any_node(self, u):
        if self.tree.is_sample(u):
            leaves = [u]
        else:
            leaves = list(self.tree.leaves(u))
        return leaves

    def _update_matrices_and_output_ibd_for_sample_pair(
            self, id1, id2, current_tmrca):

        if id1 < id2:
            row, col = id2, id1
        else:
            row, col = id1, id2

        last_tmrca = self.tmrca_matrix[row, col]
        self.tmrca_matrix[row, col] = int(current_tmrca)
        last_break_pos = self.breakpos_matrix[row, col]
        self.breakpos_matrix[row, col] = self.tree.interval[0]

        assert(last_tmrca != current_tmrca)

        IBD_seg_length = self.tree.interval[0] - last_break_pos
        self.longest = max(self.longest, IBD_seg_length)

        if IBD_seg_length > self.ibd_threshold_in_bp:
            self.id1_list.append(row)
            self.id2_list.append(col)
            self.start_list.append(int(last_break_pos))
            self.end_list.append(int(self.tree.interval[0]))
            self.length_list.append(IBD_seg_length / self.bp_per_cM)
            self.tmrca_list.append(last_tmrca)
            self.num_ibd += 1

    def output_last_ibds_for_all_pairs(self):
        for col in range(0, self.num_samples - 1):
            for row in range(col + 1, self.num_samples):

                tmrca = self.tmrca_matrix[row, col]
                last_break_pos = self.breakpos_matrix[row, col]
                IBD_seg_length = self.tree.interval[1] - last_break_pos

                if IBD_seg_length > self.ibd_threshold_in_bp:
                    self.id1_list.append(row)
                    self.id2_list.append(col)
                    self.start_list.append(int(last_break_pos))
                    self.end_list.append(int(self.tree.interval[0]))
                    self.length_list.append(
                        IBD_seg_length / self.bp_per_cM)
                    self.tmrca_list.append(tmrca)
                    self.num_ibd += 1

    def _update_matrices_due_to_new_ancestral_node(self, u):
        assert(not self.tree.is_sample(u))
        current_tmrca = self.tree.time(u)
        left_child = self.tree.left_child(u)
        right_child = self.tree.right_child(u)
        left_child_leaves = \
            self._find_tree_leaves_of_any_node(left_child)
        right_child_leaves = \
            self._find_tree_leaves_of_any_node(right_child)
        for id1 in left_child_leaves:
            for id2 in right_child_leaves:
                self._update_matrices_and_output_ibd_for_sample_pair(
                    id1, id2, current_tmrca)
        num_pairs_updated = len(left_child_leaves) * len(right_child_leaves)

        return num_pairs_updated

    def _update_matries_due_to_new_leaves(self, M, RN):

        new_leaves = self._find_tree_leaves_of_any_node(M)
        run_node = M
        num_pairs_updated = 0
        while run_node != RN:
            parent = self.tree.parent(run_node)
            current_tmrca = self.tree.time(parent)
            if run_node == self.tree.left_child(parent):
                sib = self.tree.right_child(parent)
            else:
                sib = self.tree.left_child(parent)
            other_side_leaves = self._find_tree_leaves_of_any_node(sib)
            for id1 in new_leaves:
                for id2 in other_side_leaves:
                    self._update_matrices_and_output_ibd_for_sample_pair(
                        id1, id2, current_tmrca)
            num_pairs_updated += len(new_leaves) * len(other_side_leaves)
            run_node = parent

        return num_pairs_updated

    def _update_tmrca_matrix_due_to_new_tree(self, tree):
        num_pairs_updated = 0
        for u in tree.nodes():
            # access internal nodes
            if not tree.is_sample(u):
                num_pairs_updated += \
                    self._update_matrices_due_to_new_ancestral_node(u)
        return num_pairs_updated

    def _find_edge_diff_pattern(self, diff, tree):
        """
        Most complicated algorithm here: is to find the key nodes from the
        edge_diffs record
        """

        # build parent to child and child to parent mapping for edges_out and
        # for edges_in
        edges_out = diff[1]
        edges_in = diff[2]
        in_cp_map = {}
        in_pc_map = {}
        out_cp_map = {}
        out_pc_map = {}

        for edge in edges_in:
            in_cp_map[edge.child] = edge.parent
            if edge.parent not in in_pc_map:
                in_pc_map[edge.parent] = [edge.child]
            else:
                in_pc_map[edge.parent].append(edge.child)

        for edge in edges_out:
            out_cp_map[edge.child] = edge.parent
            if edge.parent not in out_pc_map:
                out_pc_map[edge.parent] = [edge.child]
            else:
                out_pc_map[edge.parent].append(edge.child)

        num_edges = len(diff[1])

        # ========== Find the N and M node ===============
        # condition 1: two-edge pattern: elevating root,
        # can also explained as moving node breaks into edge of the root
        if num_edges == 2:  # only two edges
            N = edges_in[0].parent
            M = in_pc_map[N][0]  # M is any child of N
            RN = N
            pat = "Pattern 1: two-edge, elevating root"

        # condition 2: three-edge patterns
        elif num_edges == 3:
            two_children_Nodes = [parent for parent in in_pc_map
                                  if len(in_pc_map[parent]) == 2]
            two_children_Nodes_out = [parent for parent in out_pc_map
                                      if len(out_pc_map[parent]) == 2]
            # find N node:
            Ns = [node for node in two_children_Nodes
                  if node in in_cp_map or node == tree.root]
            assert(len(Ns) == 1)
            N = Ns[0]
            Pn = tree.parent(N)

            # this is specific for pattern 2B
            if len(two_children_Nodes_out) == 1:
                R = two_children_Nodes_out[0]
            else:
                R = None

            # # condition 2A: N node is on the edge of the prev tree root; thus
            # prev tree root is larger number of N node's children of the new
            # tree, the M node the smaller number

            if N == tree.root:
                M = min(in_pc_map[N])
                RN = N
                pat = "Pattern 2A: three-edge, moving to the root"

            # # condition 2B: N node is along the edge of R node in old tree; N
            # node in new tree  and R in old tree have same children.
            # the following logicials are 1)N is root of new tree 2) R has two
            # children in edges_out 3) R has a parent in edges_out 4) N share
            # children
            elif N != tree.root \
                    and R is not None \
                    and R in out_cp_map \
                    and out_cp_map[R] == in_cp_map[N]:
                assert(len(in_pc_map[N]) == 2)
                M = in_pc_map[N][0]  # M is any child of N
                RN = N  # no need to update for nodes older than N, thus RN=N
                pat = "Pattern 2B: three-edge, elevating along an edge"

            # # condition 2C: similar to four-edge pattern (only diff: moving
            # node is a child node of the root
            else:
                B_set = list(
                    set(out_pc_map[Pn]).intersection(set(in_pc_map[N])))
                assert(len(B_set) == 1)

                # M node is N's other child, B's sibling
                N_children = set(in_pc_map[N])
                assert(len(N_children.difference(B_set)) == 1)
                M = N_children.difference(B_set).pop()

                # R is the root of old tree, movement of M make R'sone child
                # disapear; R's other child , also M's sibling of the old tree
                # is the root of new tree. RN is root of the new tree
                RN = tree.root
                pat = "Pattern 2C: three-edge, moving a root child to an edge"

        # condition 3: four-edge pattern
        else:  #
            # N node
            two_children_Nodes = [parent for parent in in_pc_map
                                  if len(in_pc_map[parent]) == 2]
            N_list = [node for node in two_children_Nodes if node in in_cp_map]
            assert(len(N_list) == 1)
            N = N_list[0]

            # M node is N's other child, B's sibling
            Pn = tree.parent(N)
            B_set = list(set(out_pc_map[Pn]).intersection(set(in_pc_map[N])))
            assert(len(B_set) == 1)
            N_children = set(in_pc_map[N])
            assert(len(N_children.difference(B_set)) == 1)
            M = N_children.difference(B_set).pop()
            # RN node
            R = out_cp_map[M]  # R is parent of M in old tree
            Rp = out_cp_map[R]
            RN = tree.mrca(Rp, N)
            pat = "Pattern 3: four-edge, general"

        # print("RN: {}. Is root? {}".format(RN, RN == tree.root))

        return M, N, RN, pat

    def get_ibd_df(self):
        df_ibd = pd.DataFrame({'id1': self.id1_list,
                               'id2': self.id2_list,
                               'start': self.start_list,
                               'end': self.end_list,
                               'length': self.length_list,
                               'tmrca': self.tmrca_list})
        return df_ibd


if __name__ == "__main__":

    # simulation
    tree_sequence = msprime.simulate(
        sample_size=100, Ne=1000, length=30000000, recombination_rate=1e-8)

    print("Done with simulation!", file=sys.stderr)

    # meta info
    meta_info = {"chrom": 1, "bp_per_cM": 1000000,
                 "ibd_threshold_in_bp": 2000000}

    # run ts_ibd
    test_tsibd = ts_ibd(tree_sequence, meta_info)
    test_tsibd.initialize_tmrca_matrix_from_first_tree()
    test_tsibd.iterate_trees_to_update_matrixes()
    test_tsibd.output_last_ibds_for_all_pairs()
    df_ibd = test_tsibd.get_ibd_df()

    print(df_ibd)
