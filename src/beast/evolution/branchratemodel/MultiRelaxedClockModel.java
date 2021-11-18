package beast.evolution.branchratemodel;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.MathException;

import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.MRCAPrior;
import beast.math.distributions.ParametricDistribution;

@Description("Clock model that has different strict clocks for different clades, assumes clades are monophyletic")
public class MultiRelaxedClockModel extends BranchRateModel.Base implements MultiClock, Loggable {
    public Input<ParametricDistribution> rateDistInput = new Input<ParametricDistribution>("distr", "the distribution governing the rates among branches. Must have mean of 1. The clock.rate parameter can be used to change the mean rate.", Input.Validate.REQUIRED);
    public Input<IntegerParameter> categoryInput = new Input<IntegerParameter>("rateCategories", "the rate categories associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    public Input<Tree> treeInput = new Input<Tree>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
    public Input<Boolean> normalizeInput = new Input<Boolean>("normalize", "Whether to normalize the average rate (default false).", false);

    
    final public Input<Integer> numberOfDiscreteRates = new Input<>("numberOfDiscreteRates", "the number of discrete rate categories to approximate the rate distribution by. A value <= 0 will cause the number of categories to be set equal to the number of branches in the tree. (default = -1)", -1);

    RealParameter meanRate;
    int [] map;
    boolean initialised = false;
    List<MRCAPrior> calibrations = new ArrayList<MRCAPrior>();

    int nrOfRates;
    
    @Override
    public void initAndValidate() {

        tree = treeInput.get();

        categories = categoryInput.get();
        int nCategoryCount = tree.getNodeCount() - 1;
        categories.setDimension(nCategoryCount);
        Integer[] iCategories = new Integer[nCategoryCount];
        nrOfRates = (numberOfDiscreteRates.get() > 0 ? numberOfDiscreteRates.get() : tree.getNodeCount() - 1);
        for (int i = 0; i < nCategoryCount; i++) {
            iCategories[i] = i % nrOfRates;
        }
        IntegerParameter other = new IntegerParameter(iCategories);
        categories.assignFromWithoutID(other);
        categories.setLower(0);
        categories.setUpper(nrOfRates - 1);

        distribution = rateDistInput.get();

        rates = new double[nrOfRates];
        storedRates = new double[nrOfRates];
        for (int i = 0; i < rates.length; i++) {
            try {
				rates[i] = distribution.inverseCumulativeProbability((i + 0.5) / rates.length);
			} catch (MathException e) {
				throw new IllegalArgumentException(e);
			}
        }
        System.arraycopy(rates, 0, storedRates, 0, rates.length);
        normalize = normalizeInput.get();

        meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }
        
        // pick up constraints in m_initial tree
        for (final Object plugin : tree.getOutputs()) {
            if (plugin instanceof MRCAPrior && !calibrations.contains(plugin) ) {
            	if (((MRCAPrior) plugin).isMonophyleticInput.get()) {
            		calibrations.add((MRCAPrior) plugin);
            	} else {
            		Log.warning.println("Calibration that is not monophyletic found " + ((BEASTInterface) plugin).getID());
            	}
            }
        }
        if (tree.m_initial.get() != null) {
            for (final Object plugin : tree.m_initial.get().getOutputs()) {
                if (plugin instanceof MRCAPrior && !calibrations.contains(plugin) ) {
                	if (((MRCAPrior) plugin).isMonophyleticInput.get()) {
                		calibrations.add((MRCAPrior) plugin);
                	} else {
                		Log.warning.println("Calibration that is not monophyletic found " + ((BEASTInterface) plugin).getID());
                	}
                }
            }
        }
        
        meanRate.setDimension(calibrations.size() + 1);
        
        for (int i = 0; i < calibrations.size(); i++) {
        	Log.info.println(meanRate.getID() + (i+1) + " = " + calibrations.get(i).getID() + " rate");
        }
    	Log.info.println(meanRate.getID() + meanRate.getDimension() + " = root rate");
        

        try {
            double mean = rateDistInput.get().getMean();
            if (Math.abs(mean - 1.0) > 1e-6) {
                Log.warning.println("WARNING: mean of distribution for relaxed clock model is not 1.0.");
            }
        } catch (RuntimeException e) {
            // ignore
        }
//        initialise = initialiseInput.get();
    }

    public double getRateForBranch(Node node) {
		if (!initialised) {
			map = initialise(calibrations);
			initialised = true;
		}
        if (node.isRoot()) {
            // root has no rate
            return 1;
        }
        if (recompute) {
            prepare();
            recompute = false;
        }
        if (renormalize) {
            if (normalize) {
                computeFactor();
            }
            renormalize = false;
        }

        int nodeNumber = node.getNr();

        if (nodeNumber == categories.getDimension()) {
            // root node has nr less than #categories, so use that nr
            nodeNumber = node.getTree().getRoot().getNr();
        }

        int rateCategory = categories.getValue(nodeNumber);

        int rateNr = (map[node.getNr()] >= 0 ? map[node.getNr()] : calibrations.size());
        
        return rates[rateCategory] * scaleFactor * meanRate.getValue(rateNr);
    }

    // compute scale factor

    private void computeFactor() {

        //scale mean rate to 1.0 or separate parameter

        double treeRate = 0.0;
        double treeTime = 0.0;

        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            if (!node.isRoot()) {
                int nodeNumber = node.getNr();
                if (nodeNumber == categories.getDimension()) {
                    // root node has nr less than #categories, so use that nr
                    nodeNumber = node.getTree().getRoot().getNr();
                }
                int rateCategory = categories.getValue(nodeNumber);
                treeRate += rates[rateCategory] * node.getLength();
                treeTime += node.getLength();

                //System.out.println("rates and time\t" + rates[rateCategory] + "\t" + node.getLength());
            }
        }
        //treeRate /= treeTime;

        scaleFactor = 1.0 / (treeRate / treeTime);


        //System.out.println("scaleFactor\t\t\t\t\t" + scaleFactor);
    }


    private void prepare() {
//    	if (initialise) {
//    		initialise();
//    		initialise = false;
//    	}
        //System.out.println("prepare");

        categories = (IntegerParameter) categoryInput.get();

        distribution = rateDistInput.get();

        tree = treeInput.get();

        rates = new double[nrOfRates];
        try {
            for (int i = 0; i < rates.length; i++) {
                rates[i] = distribution.inverseCumulativeProbability((i + 0.5) / rates.length);
            }
        } catch (Exception e) {
            // Exception due to distribution not having  inverseCumulativeProbability implemented.
            // This should already been caught at initAndValidate()
            e.printStackTrace();
            System.exit(0);
        }

        //if (normalize) computeFactor();
    }

    @Override
    protected boolean requiresRecalculation() {
        recompute = false;
        renormalize = true;

//        if (treeInput.get().somethingIsDirty()) {
//        	recompute = true;
//            return true;
//        }
        // rateDistInput cannot be dirty?!?
        if (rateDistInput.get().isDirtyCalculation()) {
            recompute = true;
            return true;
        }
        // NOT processed as trait on the tree, so DO mark as dirty
        if (categoryInput.get().somethingIsDirty()) {
            //recompute = true;
            return true;
        }
        if (meanRate.somethingIsDirty()) {
            return true;
        }

        return recompute;
    }

    @Override
    public void store() {
        System.arraycopy(rates, 0, storedRates, 0, rates.length);
        storedScaleFactor = scaleFactor;
        super.store();
    }

    @Override
    public void restore() {
        double[] tmp = rates;
        rates = storedRates;
        storedRates = tmp;
        scaleFactor = storedScaleFactor;
        super.restore();
    }

    ParametricDistribution distribution;
    IntegerParameter categories;
    Tree tree;

    private boolean normalize = false;
    private boolean recompute = true;
    private boolean renormalize = true;

    private double[] rates;
    private double[] storedRates;
    private double scaleFactor = 1.0;
    private double storedScaleFactor = 1.0;

	@Override
	public void init(PrintStream out) {
		String id = meanRate.getID();
        for (int i = 0; i < calibrations.size(); i++) {
        	out.append(id + "." + calibrations.get(i).getID().replaceAll(".prior", "") + "\t");
        }
        out.append(id + ".root");
	}

	@Override
	public void log(long sample, PrintStream out) {
        for (int i = 0; i < calibrations.size(); i++) {
        	out.append(meanRate.getValue(i) + "\t");
        }		
        out.append(meanRate.getValue(calibrations.size()) + "\t");
	}

	@Override
	public void close(PrintStream out) {
	}

}
