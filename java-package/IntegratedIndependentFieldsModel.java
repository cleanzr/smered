import java.util.Vector;
import java.util.Random;
import cc.mallet.types.Dirichlet;
import cc.mallet.util.Randoms;
import cc.mallet.util.Maths;

// model the P(x) as a bunch of draws from various multinomials (one per field)
// model the distortion P(y|x) as y = z, with some probability, else y is uniform,
// done on a per-field basis.
public class IntegratedIndependentFieldsModel {

  // probability model P(y).
  protected int[][]  p_y_;
  // distortion model P(x|y).
  protected int[]    p_xy_;
  
  protected int[][]  n_y_;
  protected int[]    n_xy_;

  protected int r_;

  protected Vector<FieldDescriptor> schema_;

  public IntegratedIndependentFieldsModel(Vector<FieldDescriptor> schema) {
    schema_ = schema;
    p_y_ = new double[schema.size()][];
    p_xy_ = new double[schema.size()];
    n_y_ = new int[schema.size()][];
    n_xy_ = new int[schema.size()];
    for(int i = 0; i < schema_.size(); i++) {
      p_y_[i] = new double[schema_.elementAt(i).numValues()];
      n_y_[i] = new int[schema_.elementAt(i).numValues()];
    }
    r_ = 0;
    regenProbs();
  }

  public double logLikelihood(Individual ind) {
    double p = 0.0;
    for(int i = 0; i < schema_.size(); i++) {
      p += Math.log(p_y_[i][ind.getIndividualRecord().getField(i).getValue()]);
    }
    for(LinkedRecord lr : ind.getLinkedRecords()) {
      for(int i = 0; i < schema_.size(); i++) {
        if(lr.getDistortion(i)) {
          p += Math.log(p_xy_[i]);
          p += Math.log(p_y_[i][lr.getRecord().getField(i).getValue()]);
        } else {
          p +=Math.log(1.0-p_xy_[i]);
        }
      }
    }
    return p;
  }

  public void regenProbs() {
    for(int i = 0; i < schema_.size(); i++) {
      double tot = 0.0;
      for(int j = 0; j < scema_.elementAt(i).numValues(); j++) {
        tot += 1.0 + n_y_[i][j];
      }
      for(int j = 0; j < scema_.elementAt(i).numValues(); j++) {
        p_y_[i][j] = (1.0 + n_y[i][j]) / tot;
      }
      p_x_[j] = (n_xy_[j] + 0.01) / (r_ + 1.01);
    }
  }

  public void update(Individual ind) {
    Record r = ind.getIndividualRecord();
    for(int f = 0; f < schema_.size(); f++) {
      n_y_[f][r.getField(f).getValue()]++;
    }
    for(LinkedRecord lr : ind.getLinkedRecords()) {
      r_++;
      for(int f = 0; f < schema_.size(); f++) {
        if(lr.getDistortion(f)) {
          n_y_[f][lr.getRecord().getField(f).getValue()]++;
          n_xy_[f]++;
        }
      }
    }
    regenProbs();
  }

}
