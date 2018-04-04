// Specialized version of MHSampler_beka_matches_master.java
// forking off from master on 7.29.2013
// specialized version of error rates (of master) that attempts to compute FN/FP rates and true false links for waves (already doing 
// so for overall scenario)

//added on 7.12.13 - then became version 4
  //capability to handle k-way files 
  //search for the comment 'CHANGED_kfiles' and 'PRE_kfiles'
  //must first add a folder named 'records' inside your 'data' folder
  //this 'records' folder should contain the data files you want to run the sampler on (and nothing else)

//added on 7.12.13 - then became version 5
  //capability to ask the user to input the number of innerIterations, thinIterations,
   //burnIn and maxOuterIterations.
  //search for the comment 'CHANGED_inputiterations' and 'PRE_inputiterations'

//added on 7.12.13 - then became version 6
  //capability to time the sampler
  //search for the comment 'CHANGED_timer'

//version 2 creates posterior matching function for two-way matches
//and three-way matches (which preserves transitivity); output to .csv files
//DONE: write function for two-way matches
//DONE: write function for three-way matches
//DONE: add in capability to handle k-way files (post version 2)

//import java.util.Scanner;       //CHANGED_inputiterations
import java.io.*;
import java.util.*;

public class MHSampler_beka2_error_rates_v2 {
    protected IndependentFieldsModel       model_;
    protected Random                       rand_;   
    protected Vector<Vector<LinkedRecord>> blocks_;
    protected Vector<Vector<LinkedRecord>> linked_files_;
    // keys_ contains the ground-truth unique IDs
    protected Vector<Vector<String>>       keys_;
    // keys_counts_ counts how often each unique ID actually occurs ---
    // useful for determining how many links there ought to be
    protected HashMap<String,Integer>      keys_counts_;
    // keys_matches_ maps from unique IDs to records which actually
    // use that ID
    protected HashMap<String,HashSet<LinkedRecord>> keys_matches_;
    protected Individual                   list_front_;
    protected int                          next_id_;
    protected int[]                        pattern_counts_;
    // true_links_ = number of claimed links between records
    // which match ground truth, i.e., true positives
    protected Long                          true_links_ = 0L;
    // false_links_ = number of claimed links between records
    // which do NOT match ground truth, i.e., false positives
    protected Long                          false_links_ = 0L;
    // missings_links_ = number of actual (ground-truth) links between records
    // which are not reported, i.e., false negatives
    // Not counting true negatives since (a) there are lots of them and (b)
    // can find their number from true_links_, false_links_ and missing_links_
    protected Long                          missing_links_ = 0L;
    // Number of actual (ground-truth) links and non-links, useful for
    // calculating error rates
    // Use Long here rather than int because we are dealing with combinatorially-large
    // numbers
    protected Long                          actual_links_ = 0L; 
    protected Long                          actual_nonlinks_ = 0L;
    // Nmax is the total number of records, summed across all files
    protected Long                         Nmax = 0L;
    // Setting missing_links_, actual_links_, actual_nonlinks_, Nmax to zero here is just a coding
    // convenience
    // All four will be set to appropriate values during initialization
    // confusion_matrix_ stores how often records which actually have a given
    // linkage pattern are classified by the model as really having another
    // linkage pattern
    // This is a 2D array, with rows and columns both labeled by patterns, and
    // with numbers (counts) in the array --- using a 2D array of ints, exploiting
    // the fact that patterns are represented in the code as ints
    protected int[][]                      confusion_matrix_;
    protected int                          largest_file_size_ = 0;
    protected int[][]                      cluster_idx_;
    protected int			   updated_file;
    protected int			   file_new = 0;
    protected int                          records_size;   //CHANGED_kfiles

  public MHSampler_beka2_error_rates_v2(Vector<Vector<Record>> records,
                   Vector<Vector<String>> keys,
                   IndependentFieldsModel model,
                   Vector<Vector<int[]>> blocks) {
    records_size = records.size();                      //CHANGED_kfiles
    pattern_counts_ = new int[1 << records_size];       //CHANGED_kfiles
    confusion_matrix_ = new int[1 << records_size][1 << records_size];
    model_ = model;
    keys_ = keys;
    // Initialize our hash maps
    keys_counts_ = new HashMap<String,Integer>();
    keys_matches_ = new HashMap<String,HashSet<LinkedRecord>>();
    // set up files.
    Individual last = null;
    next_id_ = 0;
    linked_files_ = new Vector<Vector<LinkedRecord>>();
    for(int f = 0; f < records_size; f++) {             //CHANGED_kfiles
      Vector<Record> rv = records.elementAt(f);
      Vector<LinkedRecord> lrv = new Vector<LinkedRecord>();
      if(rv.size() > largest_file_size_)
        largest_file_size_ = rv.size();
      for(int r = 0; r < rv.size(); r++) {
	  // Increment the global count of records
	  Nmax++;
        // initialize to one individual per record.
        Record rec = rv.elementAt(r);
        Individual ind = new Individual(new Record(rec), next_id_++);
        LinkedRecord lr = new LinkedRecord(ind, rec, f, r);
        lrv.add(lr);
        ind.addRecord(lr);
	// Add the current key to the HashMap of keys encountered
	String key = keys_.elementAt(f).elementAt(r); 
        // if we have seen this key before, increment the count of how often
	// we've seen it by one
	// otherwise, add the key to the HashMap with a value of 1
	if (keys_counts_.isEmpty() || !keys_counts_.containsKey(key)) {
            keys_counts_.put(key, 1);
	} else {
	    int current_count = keys_counts_.get(key);
	    keys_counts_.put(key, current_count+1);
        }

	// FLAG: Confusion Matrix Calculation
	// Update the set of records matching the present key
	// Get the set of current matches --- start one if it's empty
	// Add the record we're working with
	// Put the set back in the HashMap
	HashSet<LinkedRecord> current_matches = new HashSet<LinkedRecord>();
	if (keys_matches_.isEmpty() || !keys_matches_.containsKey(key)) {
	    current_matches.add(lr);
	    keys_matches_.put(key,current_matches);
	} else {
	    current_matches = keys_matches_.get(key);
	    current_matches.add(lr);
	    keys_matches_.put(key,current_matches);
	}
	// END FLAG

        pattern_counts_[ind.getPattern()]++;
        // set up linked list structure.
        if(last != null) {
          ind.setPrev(last);
          last.setNext(ind);
        } else {
          list_front_ = ind;
        }
        last = ind;
      }
      linked_files_.add(lrv);
    }

    // Work out how many links there actually are, from the HashMap which
    // records how often each key is used in the ground truth
    // The number of links induced by having a unique ID appear n times is
    // n(n-1)/2 = number of edges in a complete graph on n nodes
    for (String k : keys_counts_.keySet()) {
	//get the value associate with keys k, for each value calc n(n-1)/2 and 
	//add them up
	int num_key_counts = keys_counts_.get(k);
        actual_links_ += (num_key_counts)*(num_key_counts-1)/2;
    }
    // Total number of possible links is (Nmax choose 2)
    actual_nonlinks_ = Nmax*(Nmax-1)/2 - actual_links_;
    // Initially ALL links are missing
    missing_links_ = actual_links_; 

    // Initialize confusion matrix
    // Loop through individuals and calculate estimated and true patterns for
    // each linked record
    Individual jprime = list_front_;
    do {
	int estPattern = jprime.getPattern();
	for (LinkedRecord lr : jprime.getLinkedRecords()) {
	  String key = keys_.elementAt(lr.getFile()).elementAt(lr.getIndex());
	  Individual actualInd = new Individual(keys_matches_.get(key));
	  int actualPattern = actualInd.getPattern();
	  confusion_matrix_[estPattern][actualPattern]++;
	}
	jprime = jprime.getNext();
    } while (jprime != null);


    // set up blocking.
    blocks_ = new Vector<Vector<LinkedRecord>>();
    for(Vector<int[]> block : blocks) {
      Vector<LinkedRecord> linked_block = new Vector<LinkedRecord>();
      for(int[] index : block) {
        linked_block.add(linked_files_.elementAt(index[0]).elementAt(index[1]));
      }
      blocks_.add(linked_block);
    }
    rand_ = new Random(System.currentTimeMillis());
    cluster_idx_ = new int[records_size][largest_file_size_];             //CHANGED_kfiles
  }

  public void iter() {
    // choose a block
    Vector<LinkedRecord> block = blocks_.elementAt(rand_.nextInt(blocks_.size()));
    LinkedRecord rec1 = block.elementAt(rand_.nextInt(block.size()));
    LinkedRecord rec2 = block.elementAt(rand_.nextInt(block.size()));
    if(rec1.getFile() != rec2.getFile()) { // different files
      Individual ind1 = rec1.getIndividual();
      Individual ind2 = rec2.getIndividual();
      if(ind1 == ind2) {
        doSplit(ind1, rec1, rec2);
      } else if(ind1.isMergeable(ind2)) {
        doMerge(ind1, ind2);
      }
    } 
  }

	//when it splits or merges, it updates count of true or false links
  public void doSplit(Individual i, LinkedRecord lr1, LinkedRecord lr2) {
    double lp_old = model_.logLikelihood(i);
    Individual j = new Individual(null, 0);
    Individual k = new Individual(null, 0);
    // split the records of i into j and k.
    j.addRecord(lr1);
    k.addRecord(lr2);
    for(LinkedRecord lr : i.getLinkedRecords()) {
      if(lr != lr1 && lr != lr2) {
        (rand_.nextDouble() < 0.5 ? j : k).addRecord(lr);
      }
    }
    // copy one of the records as the individual.
    j.setIndividualRecord(new Record(j.getRecords().elementAt(rand_.nextInt(j.numRecords()))));
    k.setIndividualRecord(new Record(k.getRecords().elementAt(rand_.nextInt(k.numRecords()))));
    // gibbs sample some new distortion parameters. is this what rob means?? i think rob is really sampling y and z. 
    for(LinkedRecord lr : j.getLinkedRecords()) {
      model_.sampleDistortion(j.getIndividualRecord(), lr);
    }
    for(LinkedRecord lr : k.getLinkedRecords()) {
      model_.sampleDistortion(k.getIndividualRecord(), lr);
    }
    double lp_new = model_.logLikelihood(j) + model_.logLikelihood(k);
    if(Math.log(rand_.nextDouble()) < lp_new - lp_old) {
      j.setIndex(next_id_++);
      k.setIndex(next_id_++);
      // link records to j, k
      for(LinkedRecord lr : j.getLinkedRecords()) {
        lr.setIndividual(j);
      }
      for(LinkedRecord lr : k.getLinkedRecords()) {
        lr.setIndividual(k);
      }
      // put j, k into the list and remove i.
      j.setNext(k);
      k.setPrev(j);
      k.setNext(list_front_);
      list_front_.setPrev(k);
      list_front_ = j;
      removeIndividual(i);
      pattern_counts_[i.getPattern()]--;
      pattern_counts_[j.getPattern()]++;
      pattern_counts_[k.getPattern()]++;
      // update false/true link counts
      updateStats(i, true);
      updateStats(j, false);
      updateStats(k, false);
    } else {
      // undo everything. 
      // just re-sample the distortion variables.
      for(LinkedRecord lr : i.getLinkedRecords()) {
        model_.sampleDistortion(i.getIndividualRecord(), lr);
      }
    }
  }

  protected void removeIndividual(Individual i) {
    if(list_front_ == i) {
      list_front_ = list_front_.getNext();
      list_front_.setPrev(null);
    } else {
      i.getPrev().setNext(i.getNext());
      if(i.getNext() != null) {
        i.getNext().setPrev(i.getPrev());
      }
    }
  }

    

    // Update the counts of true, false and missing links
  public void updateStats(Individual i, boolean remove) {
      int false_link_increment = 0;
      int true_link_increment = 0;
      int missing_link_increment = 0;

      // This explains the logic of missing_link_increment
      // How does the number of missing links change by adding or removing an
      // individual?
      // If we remove an individual and their links, we may now miss some more
      // real links
      // If we add an individual and their links, we may now find what had been
      // missed real links
      // Do we change whether we miss a link between 2 records, only 1 of which
      // is part of the current individual?
      // Claim: no; removing an individual only makes us miss links between
      // the records the individual contains; adding an individual only shows
      // links between the individuals it contains

      // This explains the logic of incrementing each part of the confusion matrix
      // First, we need to get the linkage pattern for individual i
      // Every record linked (in the model) to i will have this estimated pattern
      // Next, for each linked record, we need to find its true linkage pattern
      // Then we increment the corresponding cell in the table if adding an individual,
      // decrement the cell if removing an individual
      // Is it worth pre-calculating actual linkage patterns, since they do not change?

      // keys_ is the sampler's vector of unique IDs (keys) obtained from
      // the data files. 
      // This will loop over the records linked to the
      // given individual i, and compare their unique IDs to the unique
      // IDs of those actually linked to i; if there is a match, the
      // count of true links is incremented, otherwise the count of false
      // links is incremented.
      // Note that keys_ is NOT the same as keys; the former is global to the
      // sampler, the latter is local to this function

      int estPattern = i.getPattern();

      Vector<String> keys = new Vector<String>();
      for(LinkedRecord lr : i.getLinkedRecords()) {
	  String key = keys_.elementAt(lr.getFile()).elementAt(lr.getIndex());
	  for(String k : keys) {
	      if(k.equals(key)) {
		  true_link_increment++;
		  //if k == keys, then we have caught a missing link
		  missing_link_increment--;
	      }
	      else
		  false_link_increment++;
	  }
	  keys.add(key);
	  // growing the vector of keys to compare to avoids double-counting

	  // Figure out the true linkage pattern of this linked record
	  // QUERY: Make a new method of the LinkedRecord class?
	  // Temporarily create an Individual containing all and only the records actually
	  // linked to lr, i.e., have same unique ID (key) as lr
	  // FLAG: Confusion Matrix Calculation
	  Individual actualInd = new Individual(keys_matches_.get(key));
	  int actualPattern = actualInd.getPattern();
	      
	  int current_count = confusion_matrix_[estPattern][actualPattern];
	  if (remove) {
	      confusion_matrix_[estPattern][actualPattern] = current_count-1;
	  } else {
	      confusion_matrix_[estPattern][actualPattern] = current_count+1;
	  }
	  // END FLAG

      }
      // remove asks if we remove the indivdiual in question (default is to add
      // individual 
      // Flip signs of link counts if individual i is being removed
      if(remove) {
	  false_link_increment = -false_link_increment;
	  true_link_increment = -true_link_increment;
	  missing_link_increment = -missing_link_increment;
      }
      // Update global totals of true and false links
      true_links_ += true_link_increment;
      false_links_ += false_link_increment;
      missing_links_ += missing_link_increment;
  }

	//calls updateStats
  public void doMerge(Individual i, Individual j) {
    double lp_old = model_.logLikelihood(i) + model_.logLikelihood(j);
    Individual k = new Individual(null, 0);
    // stuff all records into k.
    for(LinkedRecord lr : i.getLinkedRecords()) {
      k.addRecord(lr);
    }
    for(LinkedRecord lr : j.getLinkedRecords()) {
      k.addRecord(lr);
    }
    // copy one of the records as the individual.
    k.setIndividualRecord(new Record(k.getRecords().elementAt(rand_.nextInt(k.numRecords()))));
    // gibbs sample some new distortion parameters.
    for(LinkedRecord lr : k.getLinkedRecords()) {
      model_.sampleDistortion(k.getIndividualRecord(), lr);
    }
    double lp_new = model_.logLikelihood(k);
    double lp_back = -(k.numRecords()-2)*Math.log(2.0);
    if(Math.log(rand_.nextDouble()) < lp_new + lp_back - lp_old) {
      k.setIndex(next_id_++);
      // link records to k
      for(LinkedRecord lr : k.getLinkedRecords()) {
        lr.setIndividual(k);
      }
      // put k into the list and remove i and j.
      k.setNext(list_front_);
      list_front_.setPrev(k);
      list_front_ = k;
      removeIndividual(i);
      removeIndividual(j);
      pattern_counts_[i.getPattern()]--;
      pattern_counts_[j.getPattern()]--;
      pattern_counts_[k.getPattern()]++;
      // update false/true link counts
      updateStats(i, true);
      updateStats(j, true);
      updateStats(k, false);
    } else {
      // undo everything. 
      // just re-sample the distortion variables.
      for(LinkedRecord lr : i.getLinkedRecords()) {
        model_.sampleDistortion(i.getIndividualRecord(), lr);
      }
      for(LinkedRecord lr : j.getLinkedRecords()) {
        model_.sampleDistortion(j.getIndividualRecord(), lr);
      }
    }
  }
    

	//TODO: need to change code to make sure that transposes of sets of records are treated the same. 
    public void writeAllMatches(FileWriter fw){
	 Individual in = list_front_;
	 try{
	     do {
		 // j is an object of type LinkedRecord
		 //getLinkedRecords returns vector, and sort applies to lists
		 //List<LinkedRecord> getNewLinkedRecord = new Vector<LinkedRecords>(); 
		 Vector<LinkedRecord> theLinkedRecords = in.getLinkedRecords();
		 Collections.sort(theLinkedRecords);
		 //
		 for(LinkedRecord j : theLinkedRecords) {
		     //getLinkedRecord returns all the record objects LinkedRecord(class) j(object), 
		     //each containing file and index (among other things)
		     // print the ID number of the linked record ID j using pw.println()
		     //followed by a comma (Mike Comment: (Comma = end of a row)
		     // j.foo returns the attribute named foo of the object j 
		     //output now is csv file with each matching set separated by \n
		     file_new = j.file_ + 1;
		     fw.append(file_new + "." + j.index_ + ";");
		 }
		 fw.append("\n");
		 in = in.getNext(); // Move on to next latent individual
	     } while (in != null);  // do/while loop iterates over latent individuals
	     fw.append("\n");
	 }
 	catch(Exception e) {
	    e.printStackTrace();
	}
    }
    

  public void writeLinkage(String fp) {
    for(int i = 0; i < records_size*largest_file_size_; i++) {              //CHANGED_kfiles
      cluster_idx_[i/largest_file_size_][i%largest_file_size_] = -1;
    }
    Individual in = list_front_;
    int j = 0;
    do {
      Vector<LinkedRecord> lv = in.getLinkedRecords();
      for(LinkedRecord l : lv) {
        cluster_idx_[l.getFile()][l.getIndex()] = j;
      }
      j++;
      in = in.getNext();
    } while (in != null);
    try {
      java.io.PrintWriter pw = new java.io.PrintWriter(fp);
      for(int i = 0; i < largest_file_size_; i++) {
        String cluster_idx_string = "";                                     //CHANGED_kfiles (these 6 rows)
        for (int rs = 0; rs < records_size-1; rs++) {
          cluster_idx_string = cluster_idx_string + cluster_idx_[rs][i] + ",";
        }
        cluster_idx_string = cluster_idx_string + cluster_idx_[records_size-1][i];
        pw.println(cluster_idx_string);

        // pw.println(cluster_idx_[0][i] + "," + cluster_idx_[1][i] + "," + cluster_idx_[2][i]);      //PRE_kfiles
      }
      pw.close();
    } catch(Exception e) {
      e.printStackTrace();
    }
  } 


    public static String printPattern(int pattern, int numFiles) {
	String printable = Integer.toBinaryString(pattern);
	while (printable.length() < numFiles) { printable = "0" + printable; }
	return(printable);
    }


    //the false and true links need to be averaged (take out a burn-in)  
  public static void main(String[] argv) {
      long startTime = System.currentTimeMillis();    //CHANGED_timer
      String path = "data/records";                   //CHANGED_kfiles (these 7 lines)
      File folder = new File(path);
      String[] files = folder.list();

      for (int i = 0; i < files.length; i++) {
        files[i] = path + "/" + files[i];
      }

      // String[] files = {"data/nltcs_82.txt",           //PRE_kfiles (these 3 lines)
      // "data/nltcs_89.txt",
      // "data/nltcs_94.txt"};

      String schema_fp = "data/nltcs_schema.txt"; //ATTN
      //String[] blocking_fields = {"SEX", "DOB_YEAR"};
      String[] blocking_fields = {"SEX"}; //ATTN
      RecordLoader rl = new RecordLoader(schema_fp, files);
      RecordBlocking rb = new RecordBlocking(rl.getSchema(), blocking_fields);
      Vector<Vector<int[]>> blocking = rb.constructBlocking(rl.getRecords());
      IndependentFieldsModel pr = new IndependentFieldsModel(rl.getSchema());
      MHSampler_beka2_error_rates_v2 mh = new MHSampler_beka2_error_rates_v2(rl.getRecords(), rl.getKeys(), pr, blocking);
      pr.sampleParameters(mh.list_front_);
      pr.setDistortionParam(rl.fieldIndex("DOB_DAY"), 0.000001);
      pr.setDistortionParam(rl.fieldIndex("DOB_MONTH"), 0.000001);
      pr.setDistortionParam(rl.fieldIndex("REGOFF"), 0.000001);
      pr.setDistortionParam(rl.fieldIndex("STATE"), 0.000001); //ATTN
      
      //PRE_inputiterations (next 4 lines)
      int innerIterations = 100;   // split-merge (MH) steps per outer iteration
      int thinIterations = 1;       // Write output every so many Gibbs iterations
      int burnIn = 1;              // Begin taking averages only after this many Gibbs iterations
      int maxOuterIterations = 100 + burnIn; // Number of Gibbs iterations
      
      //STARTS SCANNER    //CHANGED_inputiterations (these 17 lines)
      //PROCESS: The following chunk of code allows your console to ask you to input the iterations, thin, and burnIn
        //creates a new Scanner named 'reader' that takes in System.in
        //System.in allows you to input from your keyboard
        //System.err.println() disguises the string as an error shown on the screen
        //nextInt() scans for your input and sets it as an integer
        //then your variable (ie. innerIterations) is set to this integer
      //Scanner reader = new Scanner(System.in);
     // System.err.println("Enter the number of innerIterations");
       // int innerIterations = reader.nextInt();  // split-merge (MH) steps per outer iterations
     // System.err.println("Enter the number of thinIterations");  
       // int thinIterations = reader.nextInt();       // Write output every so many Gibbs iterations
    //  System.err.println("Enter the number of burnIn");
        //int burnIn = reader.nextInt();             // Begin taking averages only after this many Gibbs iterations
     // System.err.println("Enter the number of maxOuterIterations (excluding the burnIn)");  
      //  int maxOuterIterations = reader.nextInt() + burnIn; // Number of Gibbs iterations
      //ENDS SCANNER 
	  //true postives are true links and false negatives (falsely identified as not being linked -- missing link)
	  //true negative == truly reported negative estimate = reported as correctly as not being linked (not counting these)
	  //true postives and false negatives = denominator

      double avgTrueLinks = 0.0;       // Time-averaged numbers of true and false record linkages
      double avgFalseLinks = 0.0;
      double avgMissingLinks = 0.0;
      double FNR = 0.0;                // False Negative Rate and its time-averaged version
      double avgFNR = 0.0;
      double FPR = 0.0;                // False Positive rate and its time-averaged version
      double avgFPR = 0.0;
      double[] avgPatternCounts; // Time-averaged vector of pattern counts
      double[][] avg_confusion_matrix;
      double t = 0.0; // Counter for number of iterations post burn-in
      // Made a floating-point number so we can divide properly
      avgPatternCounts = new double[1<<files.length];                         //CHANGED_kfiles
      for (int s=0; s < avgPatternCounts.length; s++) {                       //CHANGED_kfiles
        avgPatternCounts[s] = 0.0;
      }
      avg_confusion_matrix = new double[1 << files.length][1 << files.length];
      for (int s=0; s < avg_confusion_matrix.length; s++) {
 	for (int r=0; r < avg_confusion_matrix[s].length; r++) {
	    avg_confusion_matrix[s][r] = 0.0;
	}
      }


      // Create the uniform format for file names
      int numDigits = (int)Math.ceil(Math.log10((double)maxOuterIterations));
      String fileFormat = "samples/linkage_iter_%0" + String.format("%d",numDigits)+"d.csv";
      String fileFormat2 = "samples/matches_iter_%0" + String.format("%d",numDigits)+"d.csv";
      String fileAllMatches = "samples/matches.csv";
      // we need to call mh.writematches and make sure we grab everything

      System.out.println("Nmax: "+ mh.Nmax);
      System.out.println("Actual links: " + mh.actual_links_);
      System.out.println("Actual nonlinks: " + mh.actual_nonlinks_);
	  
      File allmatchFile = new File(fileAllMatches);
      try {
	  FileWriter fw = new FileWriter(allmatchFile);
      
	  for(int i = 0; i < maxOuterIterations; i++) {
	      // do some split/merge
	      if((i >= burnIn) && (i % thinIterations == 0)) {
			  mh.writeLinkage(String.format(fileFormat, i-burnIn));
			  //call writeAllMatches which should produce one csv file (columns are Gibbs iterations)
			  mh.writeAllMatches(fw);

	      }
	      for(int j = 0; j < innerIterations; j++) {
		  mh.iter();
	      }
	      // resample model parameters
	      pr.sampleParameters(mh.list_front_);
	      // fix certain distortion model parameters
	      pr.setDistortionParam(rl.fieldIndex("DOB_DAY"), 0.000001);
	      pr.setDistortionParam(rl.fieldIndex("DOB_MONTH"), 0.000001);
	      pr.setDistortionParam(rl.fieldIndex("REGOFF"), 0.000001);
	      pr.setDistortionParam(rl.fieldIndex("STATE"), 0.000001); //ATTN

	      // DISABLING printing out status
	      // // print out status
	      // System.out.println("---------");
	      // for(int s = 1; s < 8; s++) {
	      //  String str = Integer.toBinaryString(s);
	      //  while(str.length() < 3) { str = "0" + str; }
	      //  System.out.println(str + " " + mh.pattern_counts_[s]);
	      // }
	      // System.out.println("---------");
	
	      // This prints true and false links at each iteration
	      // We want the average, calculated after a burn-in which we will specify
	      if (i > burnIn) {
		  t = (double) i - burnIn;
			  
				//using sam's paper; divide by true pos. + false negatives = true links + missing links
		  FNR =(double) mh.missing_links_/(mh.actual_links_);
		  FPR =(double) mh.false_links_/(mh.actual_links_);

		  avgTrueLinks = avgTrueLinks*((t-1)/t) + mh.true_links_*(1/t);
		  avgFalseLinks = avgFalseLinks*((t-1)/t) + mh.false_links_*(1/t);
		  avgMissingLinks = avgMissingLinks*((t-1)/t) + mh.missing_links_*(1/t);
		  avgFNR = avgFNR*((t-1)/t) + FNR*(1/t);
		  avgFPR = avgFPR*((t-1)/t) + FPR*(1/t);
		  for (int s = 1; s < avgPatternCounts.length; s++) {                                     //CHANGED_kfiles
          avgPatternCounts[s] = avgPatternCounts[s]*((t-1)/t) + mh.pattern_counts_[s]*(1/t);
      }
		  for (int s= 1; s< avg_confusion_matrix.length; s++) {
		      for (int r=1; r<avg_confusion_matrix[s].length; r++) {
			  avg_confusion_matrix[s][r] = avg_confusion_matrix[s][r]*((t-1)/t) + mh.confusion_matrix_[s][r]*(1/t);
		      }
		  }
		  System.out.println("AVG TRUE LINKS: " + avgTrueLinks);
		  System.out.println("AVG FALSE LINKS: " + avgFalseLinks);
		  System.out.println("AVG MISSING LINKS: " + avgMissingLinks);
		  System.out.println("AVG FNR: " + avgFNR);
		  System.out.println("AVG FPR: " + avgFPR);
		  
	      }
	      System.out.println("TRUE LINKS: " + mh.true_links_);
	      System.out.println("FALSE LINKS: " + mh.false_links_);
	      System.out.println("MISSING LINKS: " + mh.missing_links_);
	      System.out.println("FNR: " + FNR);
	      System.out.println("FPR: " + FPR);
	  }
	  try{
	      fw.close();
	  }
	  catch(Exception e) {
	      e.printStackTrace();
	  }
      } catch (Exception e) {
	  e.printStackTrace();
      }
      // Print out averages at the very end
      System.out.println("##### RESULTS! ######");
      System.out.println("TRUE LINKS: " + avgTrueLinks);
      System.out.println("FALSE LINKS: " + avgFalseLinks);
      // Print out the time-averaged pattern counts
      for (int s=1; s< avgPatternCounts.length; s++) {                //CHANGED_kfiles
	  //String str = Integer.toBinaryString(s);
          //while(str.length() < files.length) { str = "0" + str; }     //CHANGED_kfiles
	  System.out.println(printPattern(s,files.length) + " " + avgPatternCounts[s]);
      }
      // Print out the confusion matrix
      // header for confusion matrix --- actual linkage patterns
      for (int s=0; s<avg_confusion_matrix.length; s++) {
	  System.out.print(printPattern(s,files.length) + " ");
      }
      System.out.println("\n");
      // Each row of the confusion matrix begins with the estimated linkage pattern
      for (int s=1; s<avg_confusion_matrix.length; s++) {
	  System.out.print(printPattern(s,files.length) + " ");
	  // Then the actual counts, in the same order
	  for (int r=1; r<avg_confusion_matrix[s].length; r++) {
	      System.out.print(avg_confusion_matrix[s][r] + " ");
	  }
	  System.out.println("\n");
      }


      long endTime = System.currentTimeMillis();                   //CHANGED_timer (these 4 lines)
      System.out.println("##### FINISH ######");
      System.out.print("It took " + (endTime - startTime) + " milliseconds.");
      // if((endTime - startTime) > 60000){                        //Still needs to be fixed
      //   System.out.print(" Equivalent to about " + ((endTime - startTime)/3600000) + " hours and ");
      //   if(((endTime - startTime)%3600000) == 0) {
      //     System.out.println("0 minutes.");
      //   }
      //   else{
      //     System.out.println(((endTime - startTime) - (((endTime - startTime)/3600000)*360000))/60000 + " minute(s).");
      //   }
      // }


  }

  
}
