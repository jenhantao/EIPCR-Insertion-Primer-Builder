/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eipcr.insertion.primer.builder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 *
 * @author jenhan
 */
public class Controller {

    /**
     * parses a file opened up by JFileChooser for inputs
     */
    public void parseInputFile() {
        java.io.BufferedReader inFile = null;
        try {
            JFileChooser chooser = new JFileChooser();
            FileNameExtensionFilter apeFilter = new FileNameExtensionFilter("Ape File", "ape", "str");
            chooser.addChoosableFileFilter(apeFilter);
            chooser.showOpenDialog(null);
            if (chooser.getSelectedFile() != null) {
                inFile = new java.io.BufferedReader(new java.io.FileReader(chooser.getSelectedFile()));
                String topLine = inFile.readLine().trim();

                if (topLine.startsWith("LOCUS")) {
                    coordinates = new ArrayList();
                    insertionSiteNames = new ArrayList();
                    //skip all the lines above features
                    while (topLine != null && !topLine.startsWith("FEATURE")) {
                        topLine = inFile.readLine().trim();
                    }
                    topLine = inFile.readLine().trim();
                    String bottomLine = null;
                    while (topLine != null) {
//                        System.out.println("topLine: "+topLine);
                        if (!topLine.startsWith("/")) {
//                            System.out.println("atempting to store");
                            bottomLine = inFile.readLine().trim();
//                            System.out.println("bottom line: "+bottomLine);
                            if (bottomLine.indexOf("insertion_site") > 0) {
//                                System.out.println("storing");
                                coordinates.add(topLine.substring(topLine.lastIndexOf(" ") + 1));
                                insertionSiteNames.add(bottomLine.substring(bottomLine.indexOf("=") + 1));
                            }

                        }
                        topLine = inFile.readLine().trim();
                        if (topLine.startsWith("ORIGIN")) {
                            break;
                        }
                    }
                    topLine = inFile.readLine();
                    while (topLine != null) {
                        if (topLine.startsWith("//")) {
                            break;
                        }
//                        System.out.println(topLine.substring(10).trim().replaceAll(" ", ""));
                        _sequence += topLine.substring(10).trim().replaceAll(" ", "");
                        topLine = inFile.readLine();
                    }


//                    for (int i = 0; i < coordinates.size(); i++) {
//                        System.out.println(insertionSiteNames.get(i) + " " + coordinates.get(i));
//                    }
//                    System.out.println(_sequence);

                } else {
                    javax.swing.JOptionPane.showMessageDialog(null, "Selected file does not appear to be a genbank file", "", JOptionPane.ERROR_MESSAGE);
                }
            }
        } catch (IOException ex) {
        } finally {
            try {
                inFile.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }

        }
    }

    public Double calcTempBasic(String seq) {
        seq = seq.toUpperCase();
        int numA = 0;
        int numT = 0;
        int numC = 0;
        int numG = 0;
        for (int i = 0; i < seq.length(); i++) {
            if (seq.substring(i, i + 1).equals("A")) {
                numA++;
            } else if (seq.substring(i, i + 1).equals("T")) {
                numT++;
            } else if (seq.substring(i, i + 1).equals("C")) {
                numC++;
            } else if (seq.substring(i, i + 1).equals("G")) {
                numG++;
            }
        }
        Double toReturn = 0.0;
        if (seq.length() < 14) {
            toReturn = new Double((numA + numT) * 2 + (numG + numC) * 4);
        } else {
            toReturn = 64.9 + new Double(41 * (numG + numC - 16.4) / (numA + numT + numG + numC));

        }

        return Math.ceil(toReturn);
    }

    public Double calcTempNN(String seq) {
        int len = seq.length();
        double concP = 100 * java.lang.Math.pow(10, -9);
        double dH = 0;
        double dS = 0;
        double logCt = 0;
        double R = 1.987;
        double temp;
        String pair;
        seq = seq.toUpperCase();

        // Checks terminal base pairs
        char init = seq.charAt(0);
        if (init == 'G' || init == 'C') {
            dH += 0.1;
            dS += -2.8;
        } else if (init == 'A' || init == 'T') {
            dH += 2.3;
            dS += 4.1;
        }
        init = seq.charAt(len - 1);
        if (init == 'G' || init == 'C') {
            dH += 0.1;
            dS += -2.8;
        } else if (init == 'A' || init == 'T') {
            dH += 2.3;
            dS += 4.1;
        }

        // Checks nearest neighbor pairs
        for (int i = 0; i < len - 1; i++) {
            pair = seq.substring(i, i + 2);
            if (pair.equals("AA") || pair.equals("TT")) {
                dH += -7.9;
                dS += -22.2;
            } else if (pair.equals("AG") || pair.equals("CT")) {
                dH += -7.8;
                dS += -21.0;
            } else if (pair.equals("AT")) {
                dH += -7.2;
                dS += -20.4;
            } else if (pair.equals("AC") || pair.equals("GT")) {
                dH += -8.4;
                dS += -22.4;
            } else if (pair.equals("GA") || pair.equals("TC")) {
                dH += -8.2;
                dS += -22.2;
            } else if (pair.equals("GG") || pair.equals("CC")) {
                dH += -8.0;
                dS += -19.9;
            } else if (pair.equals("GC")) {
                dH += -9.8;
                dS += -24.4;
            } else if (pair.equals("TA")) {
                dH += -7.2;
                dS += -21.3;
            } else if (pair.equals("TG") || pair.equals("CA")) {
                dH += -8.5;
                dS += -22.7;
            } else if (pair.equals("CG")) {
                dH += -10.6;
                dS += -27.2;
            }
        }

        // Checks for self-complementarity
        int mid;
        if (len % 2 == 0) {
            mid = len / 2;
            if (seq.substring(0, mid).equals(revComp(seq.substring(mid, len)))) {
                dS += -1.4;
            }
        } else {
            mid = (len - 1) / 2;
            if (seq.substring(0, mid).equals(revComp(seq.substring(mid + 1, len)))) {
                dS += -1.4;
            }
        }

        // dH is in kCal, dS is in Cal - equilibrating units
        dH = dH * 1000;

        // logCt = java.lang.Math.log(1 / concP);
        logCt = java.lang.Math.log(concP);

        temp = (dH / (dS + (R * logCt))) - 273.15;

        //return temp;
        return temp;
    }

    public String revComp(String seq) {
        String toReturn = "";
        for (int i = 0; i < seq.length(); i++) {
            if (seq.substring(i, i + 1).equals("A")) {
                toReturn = "T" + toReturn;
            } else if (seq.substring(i, i + 1).equals("G")) {
                toReturn = "C" + toReturn;
            } else if (seq.substring(i, i + 1).equals("C")) {
                toReturn = "G" + toReturn;
            } else if (seq.substring(i, i + 1).equals("T")) {
                toReturn = "A" + toReturn;
            }
            if (seq.substring(i, i + 1).equals("a")) {
                toReturn = "t" + toReturn;
            } else if (seq.substring(i, i + 1).equals("g")) {
                toReturn = "c" + toReturn;
            } else if (seq.substring(i, i + 1).equals("c")) {
                toReturn = "g" + toReturn;
            } else if (seq.substring(i, i + 1).equals("t")) {
                toReturn = "a" + toReturn;
            } else {
                toReturn = "!" + toReturn;
            }
        }
        System.out.println("normal: "+seq);
        System.out.println("revcomp:"+toReturn);
        return toReturn;
    }

    public String generateSpacer(int length) {
        String toReturn = "";
        for (int i = 0; i < length; i++) {
            int rand = (int) Math.ceil(Math.random() * 4);
            if (rand == 1) {
                toReturn += "A";
            } else if (rand == 2) {
                toReturn += "T";
            } else if (rand == 3) {
                toReturn += "C";
            } else if (rand == 4) {
                toReturn += "G";
            }
        }
        return toReturn;
    }

    public String generateGlyCodon() {
        String toReturn = "";
        int rand = (int) Math.ceil(Math.random() * 4);
        if (rand == 1) {
            toReturn += "GGT";
        } else if (rand == 2) {
            toReturn += "GGC";
        } else if (rand == 3) {
            toReturn += "GGA";
        } else if (rand == 4) {
            toReturn += "GGG";
        }
        return toReturn;
    }

    public String generateSerCodon() {
        String toReturn = "";
        int rand = (int) Math.ceil(Math.random() * 4);
        if (rand == 1) {
            toReturn += "TCT";
        } else if (rand == 2) {
            toReturn += "TCC";
        } else if (rand == 3) {
            toReturn += "TCA";
        } else if (rand == 4) {
            toReturn += "TCG";
        }
        return toReturn;
    }

    public void generatePrimers() {
        int oligoNumber = 0; //counter for naming the oligos
//        int readingFrame = _sequence.toUpperCase().indexOf("ATG") % 3; //calculates the readingFrame
        ArrayList<String> output = new ArrayList();
        for (int i = 0; i < coordinates.size(); i++) {
            String currentCoordinates = coordinates.get(i);
            int start = Integer.parseInt(currentCoordinates.substring(0, currentCoordinates.indexOf(".")));
            int end = Integer.parseInt(currentCoordinates.substring(currentCoordinates.lastIndexOf(".") + 1));
            System.out.println("Looking at" + insertionSiteNames.get(i) + ", which is located at (" + start + "," + end + ")");
            System.out.println(_sequence.substring(start - 1, end));
            String currentSite = _sequence.substring(start - 1, end);
            for (int j = 1; j < currentSite.length() / 3; j++) {
                //j is the splitting point
                System.out.println(currentSite.substring(0, j * 3) + " " + currentSite.substring(j * 3));
                String forwardGS = generateGlyCodon() + generateSerCodon();
                String reverseGS = generateGlyCodon() + generateSerCodon();
                //The gly ser linker is the overhang that will be generated by BsaI;
                while (forwardGS.substring(1).equals(reverseGS.substring(1))) {
                    forwardGS = generateGlyCodon() + generateSerCodon();
                    reverseGS = generateGlyCodon() + generateSerCodon();
                }
                String forwardPrimer = generateSpacer(3) + "ggtctc" + forwardGS;
                String reversePrimer = generateSpacer(3) + "ctctgg" + reverseGS; //reverse primer coordiantes are hardcoded
                //homology region will be incrementally increased if Tm conditions are not met; forwardEnd and reverseEnd will be used to increment
                int forwardEnd = end;
                int reverseEnd = start - 1; //the reverse primer's homology region is upstream of the forward primer's
                String forwardHomology = _sequence.substring(start + 3 * j - 1, forwardEnd);
                String reverseHomology = _sequence.substring(reverseEnd, start + 3 * j - 1);
                System.out.println("forwardHomology: " + forwardHomology);
                System.out.println("reverseHomology: " + reverseHomology);
                while (calcTempBasic(forwardPrimer + forwardHomology) > 55 || calcTempBasic(forwardPrimer + forwardHomology) < 53) {
                    if (calcTempBasic(forwardPrimer + forwardHomology) > 55) {
                        forwardEnd--;
                    } else {
                        forwardEnd++;
                    }
                    forwardHomology = _sequence.substring(start + 3 * j - 1, forwardEnd);
                }
                while (calcTempBasic(reversePrimer + reverseHomology) > 55 || calcTempBasic(reversePrimer + reverseHomology) < 53) {
                    if (calcTempBasic(reversePrimer + reverseHomology) > 55) {
                        reverseEnd++;
                    } else {
                        reverseEnd--;
                    }
                    reverseHomology = _sequence.substring(reverseEnd, start + 3 * j - 1);
                }
                System.out.println("end fwdHomology: " + forwardHomology);
                System.out.println("end revHomology: " + reverseHomology);
                forwardPrimer += forwardHomology;
                reversePrimer = revComp(reversePrimer + reverseHomology);
                String nameStart = Integer.toString(oligoNumber);
                while (nameStart.length() < 3) {
                    nameStart = "0" + nameStart;
                }
                output.add(">" + _nameRoot + nameStart + "F | " + insertionSiteNames.get(i) + " | " + calcTempBasic(forwardPrimer));
                output.add(forwardPrimer);
                output.add(">" + _nameRoot + nameStart + "R | " + insertionSiteNames.get(i) + " | " + calcTempBasic(reversePrimer));
                output.add(reversePrimer);
                oligoNumber++;
            }
        }



        BufferedWriter outWriter = null;
        try {
            JFileChooser chooser = new JFileChooser();
            chooser.showSaveDialog(null);
            File outFile = chooser.getSelectedFile();
            if (outFile != null) {
                outWriter = new BufferedWriter(new java.io.FileWriter(chooser.getSelectedFile()));
                for (int i = 0; i < output.size(); i++) {
                    outWriter.write(output.get(i));
                    outWriter.newLine();
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                outWriter.close();
            } catch (IOException ex) {
            ex.printStackTrace();
            }

        }

    }

    public void setNameRoot(String s) {
        _nameRoot = s;
    }
    ArrayList<String> coordinates;
    ArrayList<String> insertionSiteNames;
    String _sequence = "";
    String _nameRoot = "";
}
