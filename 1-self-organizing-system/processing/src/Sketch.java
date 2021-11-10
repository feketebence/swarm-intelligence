import processing.core.PApplet;

import java.io.PrintWriter;

public class Sketch extends PApplet {

    Simulation simulation;
    PrintWriter statsFile;

    public void settings() {
        size(400, 300, P2D);
    }

    public void setup() {
        this.statsFile = createWriter("statistics/statistics.txt");

        background(255);
        fill(0);
        simulation = new Simulation(this);
    }

    public void draw() {
        background(255);
        simulation.update();
        simulation.display();

        saveFrame("frames/frame####.png");
        if (simulation.time == simulation.totalTime) {
            System.out.println("Writing into statistics file.");

            for (int i = 0; i < simulation.totalTime; i++) {
                statsFile.printf("%d %f %f %f\n", i, simulation.allMovementXRed[i], simulation.allMovementXBlue[i], simulation.allMovement[i]);
            }
            statsFile.flush();
            statsFile.close();
            System.out.println(simulation.time + " timesteps completed. Exiting...");
            exit();
        }

    }

}
