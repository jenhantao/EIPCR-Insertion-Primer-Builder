/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eipcr.insertion.primer.builder;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
        try {
            java.io.BufferedReader inFile = null;
            JFileChooser chooser = new JFileChooser();
            FileNameExtensionFilter genbankFilter = new FileNameExtensionFilter("GenBank File", "gen", "gb", "gbank", "genbank");
            chooser.addChoosableFileFilter(genbankFilter);
            chooser.showOpenDialog(null);
            if (chooser.getSelectedFile() != null) {
                inFile = new java.io.BufferedReader(new java.io.FileReader(chooser.getSelectedFile()));
                String line = inFile.readLine();
                if (line.startsWith("LOCUS")) {
                    ArrayList<String> featureLines = new ArrayList(); //holds lines that contain feature information; used to generate new features after sequence is parsed
                    String area = "";
                    String toComments = "";
                    String sequence = "";
                    while (line != null) {
                        if (line.startsWith("   ")) {
                            if (area.equals("COMMENT")) {
                                toComments = toComments + line.substring(12, line.length());
                            }
                        }
                        if (line.startsWith("COMMENT")) {
                            String comment = line.substring(12, line.length());
                            area = "COMMENT";
                            // FIXME
                            // Do something with ApE methylation data if present
                            if (comment.startsWith("ApEinfo:methylated")) {
                                String meth = line.substring(comment.length() - 1, comment.length());
                            } else {
                                toComments = toComments + comment;
                            }
                        }
                        if (line.startsWith("FEATURES")) {
                            area = "FEATURES";
                            line = inFile.readLine().trim();
                            while (!(line.startsWith("ORIGIN"))) {
                                if (!(line.startsWith("//")) && !(line.startsWith("SOURCE"))) {
                                    System.out.println(line);
                                    featureLines.add(line);
                                }
                                line = inFile.readLine().trim();
                            }
                        }
                        if (line.startsWith("ORIGIN")) {
                            line = inFile.readLine().trim();
                            while (!(line.startsWith("//"))) {
                                ArrayList<String> seq = new ArrayList(Arrays.asList(line.split(" ")));
                                for (int i = 1; i < seq.size(); i++) {
                                    sequence = sequence + seq.get(i);
                                }
                                line = inFile.readLine().trim();
                            }
                        }
                        line = inFile.readLine();
                    }

                } else {
                    javax.swing.JOptionPane.showMessageDialog(null, "Selected file does not appear to be a genbank file", "Sequence View Import", JOptionPane.ERROR_MESSAGE);
                }
                inFile.close();
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
}
