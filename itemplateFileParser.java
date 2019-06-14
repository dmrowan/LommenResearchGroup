/**
 * @author: Zaynab Ghazi
 * This program parses itemplate files and outputs
 * files containing specific data for MP and IP:
 * [ phase , phase-errors , FWHM ]
 */
import java.util.*;
import java.io.*;

public class itemplateFileParser{
    //declare indices of file columns as CONSTANTS
    public static final int TYPE=0;
    public static final int NUM=2;
    public static final int ERROR=4;
    /* this method parses an itemplate file and
     * generates 3 plain .gauss files containing the types, values
     * and errors from the itemplate.gauss
     * @param name of the itemplate file to be parsed 
     */
    public static void itemplateFileParser(String filename) throws FileNotFoundException{
	Scanner parser = new Scanner(new File(filename));
	ArrayList<String> types = new ArrayList<>();
	ArrayList<String> values = new ArrayList<>();
	ArrayList<String> errors = new ArrayList<>();

	//skip first two obsolete lines of header
	String obs1 = parser.nextLine();
	String obs2 = parser.nextLine();
	
	while (parser.hasNextLine()){
	    String line = parser.nextLine();
	    //exclude last line starting with "-"
	    if(!line.startsWith("-")){
		//split line using space as delimiter
		String[] fields = line.split(" ");
		types.add(fields[0]);
		values.add(fields[2]);
		errors.add(fields[4]);
	    }
	}
	//create to-be-written .gauss file containing values from itemplate
	PrintWriter out1 = new PrintWriter("values_"+filename);
	for(int i=0; i<values.size();i++){
	    out1.println(values.get(i));
	}
	out1.close();
	//create to-be-written .gauss file containing types from itemplate
	PrintWriter out2 = new PrintWriter("types_"+filename);
	for(int i=0; i<types.size();i++){
	    out2.println(types.get(i));
	}
	out2.close();
	//create to-be-written .gauss file containing errors from itemplate
	PrintWriter out3 = new PrintWriter("err_"+filename);
	for(int i=0; i<errors.size();i++){
	    out3.println(errors.get(i));
	}
	out3.close();
    }
    
    /*
     * main driver for the java script
     */
    public static void main(String[] args) throws FileNotFoundException,IOException{
	//Parse all itemplate files specified in the command-line
	for(int i=0; i< args.length;i++){
	    itemplateFileParser(args[i]);
	}
	/* ========= ========= MAIN_PULSE ========== =========
	 * The code below parses the resulting output files
	 *  containing the values and fetches phase and FWHM of main pulse
	 */
	PrintWriter mPhaseFile = new PrintWriter("mainPhase.gauss");
	PrintWriter mFWHMFile = new PrintWriter("mainFWHM.gauss");
	PrintWriter mPhaseError = new PrintWriter("mainPhase_errors.gauss");
	for(int i=0; i<args.length;i++){
	    String filename = "values_"+args[i];
	    String filename_="err_"+args[i];
	    double phase_query=findMainPhase(filename);
	    double phaseError_query=findMainPhaseErr(filename_);
	    double FWHM_query= findMainFWHM(filename);
	    //record energy range and corresponding phase
	    //the energy-range is a substring of the filename
	    mPhaseFile.println(args[i].substring(10,11)+"."+args[i].substring(11,13)+" "+args[i].substring(14,15)+"."+args[i].substring(15,17)+" "+phase_query);
	    //record phase_error
	    mPhaseError.println(phaseError_query);
	    //record energy range and corresponding FWHM of main pulse 
	    mFWHMFile.println(args[i].substring(10,13)+" "+args[i].substring(14,17)+" "+FWHM_query);
	    

	}
	//close written files
	mPhaseFile.close();
	mFWHMFile.close();
	mPhaseError.close();
	/* ========= ========= INTER_PULSE  ========= =========
	 * The code below parses the resulting output files
	 *  containing the values and fetches phase and FWHM of interpulse
	 */
	PrintWriter iPhaseFile_= new PrintWriter("interPhase.gauss");
	PrintWriter iPhaseError_= new PrintWriter("interPhase_errors.gauss");
	PrintWriter iFWHMFile = new PrintWriter("interFWHM.gauss");
	for(int i=0;i<args.length;i++){
	    String filename = "values_"+args[i];
	    String filename_="err_"+args[i];
	    double phase_query=findInterPhase(filename);
	    double phaseError_query=findInterPhaseErr(filename_);
	    double FWHM_query= findInterFWHM(filename);
	    //record energy range and corresponding phase
	    //the energy-range is a substring of the filename
	    iPhaseFile.println(args[i].substring(10,11)+"."+args[i].substring(11,13)+" "+args[i].substring(14,15)+"."+args[i].substring(15,17)+" "+phase_query);
	    //record phase_error
	    iPhaseError.println(phaseError_query);
	    //record energy range and corresponding FWHM of interpulse 
	    mFWHMFile.println(args[i].substring(10,13)+" "+args[i].substring(14,17)+" "+FWHM_query);
	}
	//close written files
	iPhaseFile_.close();
	iPhaseError.close();
	iFWHMFile.close();
	
    }
    /* 
     * @return phase of main pulse
     */
    public static double findMainPhase(String file) throws FileNotFoundException{
	Scanner fileParser= new Scanner(new File(file));
	//skip first 4 lines
	for (int i=0; i<4; i++) fileParser.nextLine();
	return Double.parseDouble(fileParser.nextLine());
    }
    /* 
     * @return phase-error of main pulse
     */
    public static double findMainPhaseErr(String file) throws FileNotFoundException{
	Scanner fileParser = new Scanner(new File(file));
	//skip first 4 lines
	for (int i=0; i<4; i++) fileParser.nextLine();
	return Double.parseDouble(fileParser.nextLine());
    }
    /* 
     * @return FWHM of main pulse
     */
    public static double findMainFWHM(String file) throws FileNotFoundException{
	Scanner fileParser = new Scanner (new File(file));
	//skip first 5 lines
	for (int i=0; i<5; i++) fileParser.nextLine();
	return Double.parseDouble(fileParser.nextLine());
    }
    /* 
     * @return phase of interpulse
     */
    public static double findInterPhase(String file) throws FileNotFoundException{
	Scanner fileParser= new Scanner(new File(file));
	//skip first 1 line
	fileParser.nextLine();
	return Double.parseDouble(fileParser.nextLine());
    }
    /* 
     * @return phase-error of interpulse
     */
    public static double findInterPhaseErr(String file) throws FileNotFoundException{
	Scanner fileParser = new Scanner(new File(file));
	//skip first 1 line
	fileParser.nextLine();
	return Double.parseDouble(fileParser.nextLine());
    }
    /* 
     * @return FWHM of interpulse
     */
    public static double findInterFWHM(String file) throws FileNotFpundException{
	Scanner fileParser = new Scanner (new File(file));
	//skip first 2 lines
	for (int i=0; i<2; i++) fileParser.nextLine();
	return Double.parseDouble(fileParser.nextLine());
    }
}

	

	
	

	    
	    
