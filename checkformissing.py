import bte 
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--diff_file', required=True, type=str,help='diff file')
parser.add_argument('-n', '--neighbor_file', required=True, type=str, help='output file from matUtils extract --closest-samples')
parser.add_argument('-t', '--tree_file', required=True, type=str, help='output file from matUtils extract --closest-samples')
parser.add_argument('-o', '--out_file', required=True, type=str, help='output file from matUtils extract --closest-samples')

#parser.add_argument('-t', '--threads', type=int, help='number of threads (where applicable)', default=1)
#parser.add_argument('-sd', '--scriptsDirectory', type=str, required=True, help='path to directory where scripts are')

args = parser.parse_args()

#before running this code, please run matUtils extract --within-distance --distance-threshold to extract nearest neighbors 

def read_diff(diff_file):
    #only storing missing data, stored as position and length of missing
    data = {}

    with open(diff_file) as df:
        '''
        //determine if there is a file error
        if (!file.is_open()) {
            throw std::runtime_error("Error opening file: " + diff_file);
        }
        '''
        sample = ''
        for line in df:
            if line[0] == '>':
                sample = line.strip()[1:]
                if sample in data:
                    raise Exception("Duplicate samples! Check sample: ", sample)
                else:
                    data[sample] = {}
            elif line[0] == '-':
                line = line.strip().split()
                start = int(line[1])
                length = int(line[2])
                #length = end-start
                data[sample][start] = length;    
    return data

def nearest_neighbors(neighbor_file):
    neighbs = {}
    with open(neighbor_file) as nf:
        for line in nf:
            line = line.strip().split()
            sample = line[0]
            if len(line) > 1:
                neighbors = line[1].split(',')
            else:
                neighbors = []
            neighbs[sample] = neighbors
            
    #for n in neighbs:
    #    print(n, neighbs[n])
    return neighbs

'''

def closest_samples_dfs(node, path_length, max_path_length, leaves):
    if path_length > max_path_length:
        return
    
    for child in node.children:
        if child.is_leaf():
            leaves.append((child, path_length + len(child.mutations)))
        else:
            closest_samples_dfs(child, path_length + len(child.mutations), max_path_length, leaves)

def get_closest_samples(T, nid, min_dist):
    closest_samples = ([], float('inf'))
    target = nid
    
    if not target:
        print(f"WARNING: Node {nid} not found in tree")
        return closest_samples
    
    parent = target.parent
    curr_target = target
    dist_to_orig_parent = 0
    go_up = True
    #min_dist = float('inf')

    while go_up and parent:
        parent_branch_length = len(parent.mutations) + dist_to_orig_parent
        children_and_distances = []
        min_of_sibling_leaves = float('inf')

        for child in parent.children:
            if child.is_leaf() and child.id != curr_target.id:
                child_branch_length = len(child.mutations)
                min_of_sibling_leaves = min(min_of_sibling_leaves, child_branch_length)

        for child in parent.children:
            if child.id == curr_target.id or child.id == target.parent.id:
                continue
            
            dist_so_far = dist_to_orig_parent + len(target.mutations) + len(child.mutations)
            
            if not child.is_leaf():
                leaves = []
                max_path = min_of_sibling_leaves if min_of_sibling_leaves != float('inf') else float('inf')
                closest_samples_dfs(child, dist_so_far, max_path + dist_so_far, leaves)
                print('LEAVES', leaves)
                children_and_distances.extend(leaves)
                print('CANDD', children_and_distances)
            else:
                children_and_distances.append((child, dist_so_far))
        
        for child, child_branch_length in children_and_distances:
            print('CBL', child_branch_length)
            if child_branch_length < parent_branch_length:
                print('is this happening?')
                go_up = False
            if child_branch_length < min_dist:
                min_dist = child_branch_length
                closest_samples = ([child.id], min_dist)
            elif child_branch_length == min_dist:
                closest_samples[0].append(child.id)

        curr_target = parent
        parent = curr_target.parent
        dist_to_orig_parent = parent_branch_length

    return closest_samples
'''

def LCA(tree, node1, node2):
    possible_lcas_order = [anc.id for anc in tree.rsearch(node1)]
    possible_lcas = set(possible_lcas_order)
    new_ancestors = set([anc.id for anc in tree.rsearch(node2)])
    if len(new_ancestors) == 0:
        print(f"WARNING: node {node2.id} not found in the tree! Ignoring for LCA calculations")
    possible_lcas = possible_lcas.intersection(new_ancestors)
    if len(possible_lcas) == 0:
        raise ValueError("ERROR: no valid LCA! Check that input nodes are found on the tree.")
    #if only one choice is left, just return that.
    if len(possible_lcas) == 1:
        return possible_lcas.pop()
    #otherwise, return the element of possible_lcas that's earliest in the order.
    for pl in possible_lcas_order:
        if pl in possible_lcas:
            return pl

def list_comp(mutations, missing_data):
    #this is a set because missing data is not combined so there may be some double comparisons
    #not the most efficient but at the end we will know which mutations we need to delete for each node (thats whats in delmuts)
    delmuts = set()

    #if there is more than one mutation in the list
    if len(mutations) > 1:
        #create a list for missing data (this is not efficient but its fast to code and i can delete it as i go which will gradually make it faster)
        miss = []
        #missing data for one or the other neighbor from compare data
        for p in missing_data:
            # add missing data (positon and length of missing) to miss. this gives us a mutatable varaible which we can make smaller as we iterate
            miss.append((p,missing_data[p]))
        #print('miss', miss)

        remove = set()
        for m in mutations:
            for p in miss:
                pstart = p[0]
                pend = pstart + p[1]
                #print('start', pstart, 'end', pend)
                if m > pend:
                    #remove.add(p)
                    miss.remove(p)
                    continue
                if m < pstart:
                    break
                if m >= pstart and m <= pend:
                    print("MATCH", p, m) 
                    delmuts.add(m)
                    #print(m, p)
            #for r in remove:
            #    print(r)
            #    miss.remove(r)
    #if there is only one mutation 
    elif len(mutations)==1:
        #iterate through missing data and compare 
        for m in missing_data:
            #start and end of missing data
            mstart = m
            mend = m+missing_data[m]
            #print('start', mstart, 'end', mend)

            #if the single mutation position is >mend we need to keep iterating through missing 
            if mutations[0] > mend:
                continue
            #if the single mutation is less than missing we can stop iterating
            if mutations[0] < mstart:
                break
            #if the single mutation postion overlaps the current missing position we need to record the mutation and end the loop (dont need to compare anymore bc we only have one mut and we already decided to delete it )
            if mutations[0] >= mstart and mutations[0] <= mend:
                #print("MATCH", m, mutations[0]) 
                delmuts.add(mutations[0])
                break
    return delmuts
            
            

def compare_data(missing_data, node, leaf1, leaf2, muts_to_delete):
    #print('COMPAREDATA')
    #print(len(node.mutations))

    #if node does not already have a set of deletions 
    if node.id not in muts_to_delete:  
        nodemuts = set()
    #node has set of deletions, we will add to it 
    else:
        nodemuts = set(muts_to_delete[node.id])
    #print((missing_data[leaf1]))
    #print((missing_data[leaf2]))
    
    #compare current node mutations (which may be internal or another leaf) to leaf1 and leaf2 missing data to see if the snp distance is being inflated
    if len(node.mutations) > 0:
        #mutations holds positions, strgs holds actual mutation (inefficient but it works)
        mutations = []
        strgs = {}

        #collect position and mutation itself
        for m in node.mutations:
            mutations.append(int(m[1:-1]))
            strgs[int(m[1:-1])] = (m[0], m[-1])
        mutations = sorted(mutations)
        #dont compare leaf1 to its own missing data
        #compare leaf1 to every other node in path between it and leaf2
        if node.id != leaf1:
            print('comp', leaf1, node.id)
            delmuts = list_comp(mutations, missing_data[leaf1])
            for d in delmuts:

                print('NODE', node.id, d)
                print('MUTATION', strgs[d][0]+str(d)+strgs[d][1])
                nodemuts.add(strgs[d][0]+str(d)+strgs[d][1])
        #dont compare leaf2 to its own missing data
        if node.id != leaf2:
            print('comp', leaf2, node.id)
            delmuts = list_comp(mutations, missing_data[leaf2])
            for d in delmuts:
                print('NODE', node.id, d)
                nodemuts.add(strgs[d][0]+str(d)+strgs[d][1])
            
    #if node.id not in muts_to_delete:
    #muts = []
    #for n in nodemuts:
    #    for m in node.mutations:
    #        if n == m[1:-1]:
    #            muts.append(m)

    if len(nodemuts) > 0:
        muts_to_delete[node.id] = sorted(nodemuts)
    #else:
    #    print('WHAT')
            #print('mutations', mutations)



    #if current_node.id != leaf1: 

def checknodes(tree, leaf1, leaf2, missing_data, muts_to_delete):
    #get lca of leaf and its neighbor
    lca = LCA(tree, leaf1, leaf2)
    #work up path from leaf1 to lca and compare each node to missing data
    current_node = tree.get_node(leaf1)

    #traverse the path from leaf1 to mrca
    while current_node.id != lca:
        #compare current node mutations to the missing data   
        compare_data(missing_data, current_node, leaf1, leaf2, muts_to_delete)
        #update current_node to parent
        current_node = current_node.parent
        #print(current_node.id)
        #when we hit the mrca, we move to the other side, the leaf side 
        #we dont need to include the leaf bc those mutations already dont exist?? not true? 
        
    current_node = tree.get_node(leaf2)
    #traverse path from leaf2 to mrca
    while current_node.id != lca:
        #compare current node to missing data
        #nodeComp(current_node, leaf, diff_data,missing_data);
        compare_data(missing_data, current_node, leaf1, leaf2, muts_to_delete)

        current_node = current_node.parent
        #print(current_node.id)
    
    #do the lca last
    lca = tree.get_node(lca)
    compare_data(missing_data, lca, leaf1, leaf2, muts_to_delete)

    #nodeComp(mrca, leaf, diff_data,missing_data);
    
    #this code assumes that all mutations entered are snps, will not work if mutations are longer than 1bp (which is currently usher's capability) 




#iterate through all leaves using neighbors data structure
# for each leaf iterate through all neighbors
#for each neighbor find the path between the original node and the neighbor 
#for each node in the path, compare the missing data from both leaves and identify nodes that should be deleted because of missing data 
def checkdistances(tree_file, neighbors, data):
    #global variable updates throughout 
    muts_to_delete = {}
    tree = bte.MATree(tree_file)
    #track comparisons to remove redundancy
    comps = {}
    #iterates through all leaves
    for n in neighbors:
        #if a leaf hasn't been used yet, add it to comps
        if n not in comps:
            comps[n] = []
        #for all neighbors of your current leaf n
        for s in neighbors[n]:
            #if s hasnt been compared to n (or vice versa)
            if s not in comps[n]:

                if s not in comps:
                    comps[s] = []
                #add n and s and n,s s,n to comps to prevent redundancy, they will all be compared during this iteration
                comps[s].append(n)    
                comps[n].append(s)
                checknodes(tree, n, s, data, muts_to_delete)
    total = 0
    for n in muts_to_delete:
        node = tree.get_node(n)
        parsimonyred = 1*len(muts_to_delete[n])
        total += parsimonyred
    print("TOTAL", total)


    '''
    EOD 8/5/2024
    '''
    with open(args.out_file, 'w') as of:
        for d in muts_to_delete:
            print(d)
            of.write(f'{d}\n')
            for p in muts_to_delete[d]:
                of.write(f'{p}\n')

    return muts_to_delete          
            #add to comps so im not redundant 




#retrieve missing data
data = read_diff(args.diff_file)
#parse neighbors from file, this contains all neighbors within specified distance of each leaf
neighbors = nearest_neighbors(args.neighbor_file)
#find which mutations should be deleted 
deletemuts = checkdistances(args.tree_file, neighbors, data)

