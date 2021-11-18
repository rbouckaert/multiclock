package beast.evolution.branchratemodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.MRCAPrior;


public interface MultiClock {

	default public int [] initialise(final List<MRCAPrior> clades) {
		// set map to all -1
		MRCAPrior p0 = clades.get(0);
		Tree tree = p0.treeInput.get();
		int nrOfNodes = tree.getNodeCount();
		int [] map = new int[nrOfNodes];
		Arrays.fill(map, -1);
		
		// go through the clades
		List<List<Integer>> cladeIDs = new ArrayList<>();
		for (MRCAPrior clade : clades) {
			cladeIDs.add(getCladeIDS(clade));
		}
		
		// process in order of clade size, the biggest first
		// this ensures nested clades will be processed correctly
		boolean [] done = new boolean[clades.size()];
		for (int i = 0; i < clades.size(); i++) {
			int maxLength = Integer.MIN_VALUE;
			int maxIndex = -1;
			for (int j = 0; j < clades.size(); j++) {
				if (!done[j] && cladeIDs.get(j).size() > maxLength) {
					maxLength = cladeIDs.get(j).size();
					maxIndex = j;
				}
			}
			for (Integer k : cladeIDs.get(maxIndex)) {
				map[k] = maxIndex;
			}
			done[maxIndex] = true;
		}
		return map;
	}
	
	/** returns list of node numbers of nodes in the clade **/
	default  List<Integer> getCladeIDS(MRCAPrior clade) {
	    // array of flags to indicate which taxa are in the set
	    boolean[] isInTaxaSet;
	    // array of indices of taxa
	    int[] taxonIndex;
	    int nrOfTaxa;
	    
        Tree tree = clade.treeInput.get();
        final List<String> sTaxaNames = new ArrayList<String>();
        for (final String sTaxon : tree.getTaxaNames()) {
            sTaxaNames.add(sTaxon);
        }
        // determine nr of taxa in taxon set
        List<String> set = null;
        if (clade.taxonsetInput.get() != null) {
            set = clade.taxonsetInput.get().asStringList();
            nrOfTaxa = set.size();
        } else {
            // assume all taxa
            nrOfTaxa = sTaxaNames.size();
        }
        
        // determine which taxa are in the set
        taxonIndex = new int[nrOfTaxa];
        isInTaxaSet = new boolean[sTaxaNames.size()];
        int k = 0;
        for (final String sTaxon : set) {
            final int iTaxon = sTaxaNames.indexOf(sTaxon);
            if (iTaxon < 0) {
                throw new RuntimeException("Cannot find taxon " + sTaxon + " in data");
            }
            if (isInTaxaSet[iTaxon]) {
                throw new RuntimeException("Taxon " + sTaxon + " is defined multiple times, while they should be unique");
            }
            isInTaxaSet[iTaxon] = true;
            taxonIndex[k++] = iTaxon;
        }
		List<Integer> list = new ArrayList<>();
        collectCladeNodes(tree.getRoot(), new int[1], list, isInTaxaSet, nrOfTaxa);		
		return list;
	}

	/** recurse through tree and collect nodes in a clade **/
	default int collectCladeNodes(final Node node, final int[] nTaxonCount, List<Integer> list, final boolean [] isInTaxaSet, final int nrOfTaxa) {
        if (node.isLeaf()) {
            nTaxonCount[0]++;
            if (isInTaxaSet[node.getNr()]) {
            	list.add(node.getNr());
                return 1;
            } else {
                return 0;
            }
        } else {
            int iTaxons = collectCladeNodes(node.getLeft(), nTaxonCount, list, isInTaxaSet, nrOfTaxa);
            final int nLeftTaxa = nTaxonCount[0];
            nTaxonCount[0] = 0;
            if (node.getRight() != null) {
                iTaxons += collectCladeNodes(node.getRight(), nTaxonCount, list, isInTaxaSet, nrOfTaxa);
                final int nRightTaxa = nTaxonCount[0];
                nTaxonCount[0] = nLeftTaxa + nRightTaxa;
                if (iTaxons == nrOfTaxa) {
                    return iTaxons + 1;
                }
                if (nTaxonCount[0] > 0 && nTaxonCount[0] < nrOfTaxa) {
                	list.add(node.getNr());
                }
            }
            return iTaxons;
        }
    }
}
