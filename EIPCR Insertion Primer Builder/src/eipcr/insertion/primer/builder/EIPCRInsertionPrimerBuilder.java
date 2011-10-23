
package eipcr.insertion.primer.builder;

/**
 *
 * @author jenhan
 */
public class EIPCRInsertionPrimerBuilder {

    /**
     * @param args no arguments, a JFileChooser will be used to select files
     */
    public static void main(String[] args) {
        Controller _controller = new Controller();
        _controller.parseInputFile();
        _controller.generatePrimers();
        System.out.println(_controller.calcTempBasic("GTATCACGAGGCAGAATTTCAG"));
    }
}
