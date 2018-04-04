import java.util.*;

// An individual is just a cluster of records from different files.
public class Individual {
  // The fields of this individual (using a record object for storage).
  protected Record                fields_;
  // The linked records.
  protected Vector<LinkedRecord>  records_;
  // The position in the big vector of individuals.
  protected int                   index_;
  // position in a big linked list of individuals
  protected Individual            prev_;
  protected Individual            next_;

  protected int                   cluster_index_;

  public Individual(Record fields, int index) {
    records_ = new Vector<LinkedRecord>();
    index_ = index;
    fields_ = fields;
    cluster_index_ = 0;
  }

  public Individual(Record fields, int index, int cluster_index) {
    records_ = new Vector<LinkedRecord>();
    index_ = index;
    fields_ = fields;
    cluster_index_ = cluster_index;
  }

    // FLAG: Confusion Matrix Calculation
    // Constructor that builds an individual from a set of
    // LinkedRecords --- note that it has no previous or next Individual
    // set, so this will not alter the main linked list of Individuals
    // (though we could set a previous and next if we wanted to add
    // this Individual into the linked list)
    public Individual(Set<LinkedRecord> some_records) {
	records_ = new Vector<LinkedRecord>();
	for (LinkedRecord lr : some_records) {
	    records_.add(lr);
	}
	index_ = 0;
	//fields_ = some_records[1].getRecord();
	// ATTN: Pick some arbitrary member of the set some_records and
	// use its fields
	cluster_index_ = 0;
    }

  public int getClusterIndex() {
    return cluster_index_;
  }

  public void setClusterIndex(int i) {
    cluster_index_ = i;
  }

  public Individual getNext() {
    return next_;
  }
  
  public void setNext(Individual next) {
    next_ = next;
  }

  public Individual getPrev() {
    return prev_;
  }
  
  public void setPrev(Individual prev) {
    prev_ = prev;
  }

  public void setIndex(int ind){
    index_ = ind;
  }

  public int getIndex() {
    return index_;
  }
  
  public void addRecord(LinkedRecord lr) {
    records_.add(lr);
  }

    // This method should be used with great caution --- not
    // normally wise to remove records from an individual
    public void removeRecord(LinkedRecord lr) {
	records_.remove(lr);
    }
  
  public Record getIndividualRecord() {
    return fields_;
  }

  public void setIndividualRecord(Record fields) {
    fields_ = fields;
  }

  public Vector<LinkedRecord> getLinkedRecords() {
    return records_;
  }

  public int numRecords() {
    return records_.size();
  }

  public Vector<Record> getRecords() {
    Vector<Record> recs = new Vector<Record>();
    for(LinkedRecord lr : records_) {
      recs.add(lr.getRecord());
    }
    return recs;
  }

    //bug in code here because of duplicates
    //what should the pattern counter be doing in this case
    //do increment if haven't encounter record from this file?
  public int getPattern() {
    int pattern = 0;
    for(LinkedRecord lr : records_) {
      pattern += 1 << lr.getFile();
    }
    return pattern;
  } 

  public boolean isMergeable(Individual ind) {
    return ((getPattern() & ind.getPattern()) == 0);
  }


}
