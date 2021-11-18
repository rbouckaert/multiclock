package beast.evolution.branchratemodel;



import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.MRCAPrior;

@Description("Clock model that has different strict clocks for different clades, assumes clades are monophyletic, rates are drawn from log-normal")
public class MultiRelaxedClockModel3 extends BranchRateModel.Base implements MultiClock {
    public Input<RealParameter> stdDevInput = new Input<>("stddev", "standard deviation for log normal distribution.", Input.Validate.REQUIRED);
    public Input<IntegerParameter> categoryInput = new Input<IntegerParameter>("rateCategories", "the rate categories associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    public Input<Tree> treeInput = new Input<Tree>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
    public Input<Boolean> normalizeInput = new Input<Boolean>("normalize", "Whether to normalize the average rate (default false).", false);

    RealParameter meanRate;
    int [] map;
    boolean initialised = false;
    List<MRCAPrior> calibrations = new ArrayList<>();

    LogNormalImpl distribution;
    IntegerParameter categories;
    Tree tree;

    private boolean normalize = false;
    private boolean recompute = true;
    private boolean renormalize = true;

    private double[] rates;
    private double[] storedRates;
    private double scaleFactor = 1.0;
    private double storedScaleFactor = 1.0;
    RealParameter stddevs;

    @Override
    public void initAndValidate() {

        tree = treeInput.get();

        distribution = new LogNormalImpl(1.0, 0.33);

        normalize = normalizeInput.get();

        meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }
        stddevs = stdDevInput.get();
        
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

        categories = categoryInput.get();
        int nCategoryCount = calibrations.size() + 1; //tree.getNodeCount() - 1;
        categories.setDimension(nCategoryCount);
        Integer[] iCategories = new Integer[nCategoryCount];
        for (int i = 0; i < nCategoryCount; i++) {
            iCategories[i] = i;
        }
        IntegerParameter other = new IntegerParameter(iCategories);
        categories.assignFromWithoutID(other);
        categories.setLower(0);
        categories.setUpper(tree.getNodeCount() - 1);
        
        rates = new double[tree.getNodeCount()];
        storedRates = new double[rates.length];
        for (int i = 0; i < rates.length; i++) {
            try {
				rates[i] = distribution.inverseCumulativeProbability((i + 0.5) / rates.length, stddevs.getValue());
			} catch (MathException e) {
				throw new IllegalArgumentException(e);
			}
        }
        for (int k = 0; k < rates.length; k++) {
        	System.arraycopy(rates, 0, storedRates, 0, rates.length);
        }
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


        int rateNr = (map[node.getNr()] >= 0 ? map[node.getNr()] : calibrations.size());
        int rateCategory = categories.getValue(rateNr);
        
        return rates[rateCategory] * scaleFactor * meanRate.getValue();
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
                int rateNr = (map[node.getNr()] >= 0 ? map[node.getNr()] : calibrations.size());
                int rateCategory = categories.getValue(rateNr);
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

        try {
	        for (int i = 0; i < rates.length; i++) {
	            rates[i] = distribution.inverseCumulativeProbability((i + 0.5) / rates.length, stddevs.getValue());
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
        if (stddevs.somethingIsDirty()) {
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


    public class LogNormalImpl implements ContinuousDistribution {
        double m_fMean;
        double m_fStdDev;
        NormalDistributionImpl m_normal = new NormalDistributionImpl(0, 1);

        public LogNormalImpl(double fMean, double fStdDev) {
            setMeanAndStdDev(fMean, fStdDev);
        }

        void setMeanAndStdDev(double fMean, double fStdDev) {
        	// mean in real space
            fMean = Math.log(fMean) - (0.5 * fStdDev * fStdDev);
            m_fMean = fMean;
            m_fStdDev = fStdDev;
            m_normal.setMean(fMean);
            m_normal.setStandardDeviation(fStdDev);
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
            return m_normal.cumulativeProbability(Math.log(x));
        }

        @Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            return cumulativeProbability(x1) - cumulativeProbability(x0);
        }

        public double inverseCumulativeProbability(double p, double stddev) throws MathException {
        	setMeanAndStdDev(1.0, stddev);
            return Math.exp(m_normal.inverseCumulativeProbability(p));
        }

        @Override
        public double inverseCumulativeProbability(double p) throws MathException {
            return Math.exp(m_normal.inverseCumulativeProbability(p));
        }
        
        @Override
        public double density(double fX) {
            if( fX <= 0 ) {
                return 0;
            }
            return m_normal.density(Math.log(fX)) / fX;
        }

        @Override
        public double logDensity(double fX) {
            if( fX <= 0 ) {
                return  Double.NEGATIVE_INFINITY;
            }
            return m_normal.logDensity(Math.log(fX)) - Math.log(fX);
        }
    } // class LogNormalImpl
}
