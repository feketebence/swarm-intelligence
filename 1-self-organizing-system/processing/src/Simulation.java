import processing.core.PApplet;

import java.util.Arrays;
import java.util.Random;

public class Simulation {
    PApplet p;
    Random random_generator;

    int nParticles;
    float rMin;

    float dt;
    float particleParticleScreeningLength;
    float particleDrivingForce;

    float verletCutoffDistance;
    float verletCutoffDistanceSquared;

    float verletIntershellSquared;

    float[] particles_x, particles_y;
    float[] particles_fx, particles_fy;
    float[] particles_dx_so_far, particles_dy_so_far;

    float[] particles_direction;
    int[] colors;

    int[] verletListI = null;
    int[] verletListJ = null;
    int nVerlet = 0;
    int nVerletMax = 0;

    int halfWidth, halfHeight;

    int time;
    int totalTime;
    int echoTime;

    float temperature;
    float additional_force;

    boolean flagToRebuildVerlet;

    float SX, SY, halfSX, halfSY;

    float[] allMovementXRed;
    float[] allMovementXBlue;
    float[] allMovement;

    Simulation(PApplet p) {
        this.p = p;
        this.random_generator = new Random();

        this.nParticles = 1500;
        this.rMin = 8.25f;

        this.totalTime = 2500;
        this.echoTime = 100;

        this.temperature = 0.0f;
        this.additional_force = 0.0f;

        this.dt = 1f;
        this.particleParticleScreeningLength = 10.0f;
        this.particleDrivingForce = 1.5f;

        System.out.println("Timestep (dt): " + dt);
        System.out.println("Screening length: " + particleParticleScreeningLength);
        System.out.println("Driving force on particles: " + particleDrivingForce);

        this.verletCutoffDistance = (float) 1.5 * this.particleParticleScreeningLength;
        this.verletCutoffDistanceSquared = (float) Math.pow(verletCutoffDistance, 2);
        this.verletIntershellSquared = verletCutoffDistance - particleParticleScreeningLength;
        this.verletIntershellSquared = verletIntershellSquared / 2;
        this.verletIntershellSquared = (float) Math.pow(verletIntershellSquared, 2);
        System.out.println("Verlet cutoff distance: " + verletCutoffDistance);
        System.out.println("Verlet cutoff distance squared: " + verletCutoffDistanceSquared);
        System.out.println("Half of Verlet intershell distance: " + Math.sqrt(verletIntershellSquared));
        System.out.println("Half of Verlet intershell distance squared: " + verletIntershellSquared);

        this.particles_x = new float[nParticles];
        this.particles_y = new float[nParticles];
        this.particles_fx = new float[nParticles];
        this.particles_fy = new float[nParticles];
        this.particles_dx_so_far = new float[nParticles];
        this.particles_dy_so_far = new float[nParticles];

        this.particles_direction = new float[nParticles];
        this.colors = new int[nParticles];

        this.halfWidth = p.width / 2;
        this.halfHeight = p.height / 2;

        initParticlesRandomly();
        initParticleColorsAndDirection();

        this.time = 0;


        this.SX = 120.0f;
        this.SY = 120.0f;
        this.halfSX = this.SX / 2;
        this.halfSY = this.SX / 2;

        this.allMovementXRed = new float[totalTime];
        this.allMovementXBlue = new float[totalTime];
        this.allMovement = new float[totalTime];

        // build verlet list for the first time
        rebuildVerletList();
    }

     private void initParticlesRandomly() {
        int i, j;
        float xTry, yTry;
        float dr;
        boolean overlap;
        int nTrials;

        int nParticlesSquared = (int)Math.pow(nParticles, 2);

        for(i = 0; i < nParticles; i++) {
            xTry = 0.0f;
            yTry = 0.0f;

            overlap = true;
            nTrials = 0;

            while (overlap && nTrials < nParticles) {
                xTry = this.random_generator.nextFloat() * p.width;
                yTry = this.random_generator.nextFloat() * p.height;

                overlap = false;

                for (j = 0; j < i; j++) {
                    dr = distanceFoldedPBC(xTry, yTry, particles_x[j], particles_y[j]);

                    if(dr < this.rMin) {
                        overlap = true;
                        nTrials++;
                        break;
                    }
                }
            }

            if (nTrials == nParticlesSquared) {
                System.out.println("Cannot place particles randomly, (system too dense). Quitting...");
                System.exit(1);
            }

            particles_x[i] = xTry;
            particles_y[i] = yTry;

            particles_fx[i] = 0.0f;
            particles_fy[i] = 0.0f;

            particles_fx[i] = 0.0f;
            particles_fy[i] = 0.0f;
        }

        System.out.println("Random arrangement of particles initialized, placed n = " + nParticles + " particles.");
    }

    private float distanceFoldedPBC(float x0, float y0, float x1, float y1) {
        float r;
        float dx, dy;

        dx = x1 - x0;
        dy = y1 - y0;

        if (dx > this.halfWidth)
            dx -= p.width;
        if (dx <= -this.halfWidth)
            dx += p.width;

        if (dy > this.halfHeight)
            dx -= p.height;
        if (dy <= -this.halfHeight)
            dx += p.height;

        double dxSquared = Math.pow(dx, 2);
        double dySquared = Math.pow(dy, 2);
        r = (float)Math.sqrt(dxSquared + dySquared);

        return r;
    }

/*
    void initCoordinatesUniform() {
        // init particle coordinates from a Uniform distribution
        for(int i = 0; i < this.nParticles; i++) {
            this.particles_x[i] = this.random_generator.nextFloat() * p.width;
            this.particles_y[i] = this.random_generator.nextFloat() * p.height;
//            System.out.println("particles_x[" + i + "]: " + this.particles_x[i] + ", particles_y[" + i + "]: " + this.particles_y[i]);
        }
    }

    void initCoordinatesGaussian() {
        int mid_width = p.width/2;
        int mid_height = p.height/2;

        // init particle coordinates from a Gaussian distribution
        for(int i = 0; i < nParticles; i++) {
            this.particles_x[i] = (float)random_generator.nextGaussian() * 100 + mid_width;
            this.particles_y[i] = (float)random_generator.nextGaussian() * 100 + mid_height;
//            System.out.println("particles_x[" + i + "]: " + this.particles_x[i] + ", particles_y[" + i + "]: " + this.particles_y[i]);
        }
    }
*/

    private void initParticleColorsAndDirection() {
        // init particle colors
        for(int i = 0; i < nParticles; i++) {
            if (p.random((float)0.0, (float)1.0) > 0.5) {
                this.colors[i] = 0;
                this.particles_direction[i] = -1;
            } else {
                this.colors[i] = 1;
                this.particles_direction[i] = 1;
            }
//            System.out.println("colors[" + i + "] = " + this.colors[i]);
        }
    }

    private void calculateStatistics() {
        int i;
        int nRed = 0;
        int nBlue = 0;
        double allMovementXBlue = 0.0;
        double allMovementXRed = 0.0;
        double allMovement = 0.0;

        for (i = 0; i < nParticles; i++) {
            if (colors[i] == 0) {  //RED
                allMovementXRed += particles_fx[i];
                allMovement -= particles_fx[i];
                nRed++;
            } else if (colors[i] == 1) {    // BLUE
                allMovementXBlue += particles_fx[i];
                allMovement += particles_fx[i];
                nBlue++;
            }
        }
        allMovementXRed = allMovementXRed / nRed;
        allMovementXBlue = allMovementXBlue / nBlue;
        allMovement = allMovement / (nRed + nBlue);

        if (this.time < this.totalTime) {
            this.allMovement[time] = (float) allMovement;
            this.allMovementXBlue[time] = (float) allMovementXBlue;
            this.allMovementXRed[time] = (float) allMovementXRed;
        }
    }


    private void runSimulationStep() {
        calculateExternalForcesOnParticles();
        calculatePairwiseForces();
        calculateThermalForces();

        calculateStatistics();

        moveParticles();

        checkVerletRebuildContidionAndSetFlag();

        if (this.flagToRebuildVerlet) {
            rebuildVerletList();
        }

        if (this.time % this.echoTime == 0) {
            System.out.println("Timestep: " + this.time + "/" + this.totalTime);
        }
    }

    private void calculateThermalForces() {
        int i;
        for (i = 0; i < nParticles; i++) {
            particles_fx[i] += temperature * (2.0f * random_generator.nextFloat() - 1.0f);
            particles_fy[i] += temperature * (2.0f * random_generator.nextFloat() - 1.0f);
        }
    }

    private void calculateExternalForcesOnParticles() {
        int i;

        for(i = 0; i < nParticles; i++) {
            particles_fx[i] += particles_direction[i] * this.particleDrivingForce;
        }
    }

    private void calculatePairwiseForces() {
        int i, j, ii;
        float r, r2, f;
        float dx, dy;

        for (ii = 0; ii < nVerlet; ii++) {
            i = verletListI[ii];
            j = verletListJ[ii];

            dx = particles_x[j] - particles_x[i];
            dy = particles_y[j] - particles_y[i];

            if (dx > this.halfWidth)
                dx -= p.width;
            if (dx <= -this.halfWidth)
                dx += p.width;

            if (dy > this.halfHeight)
                dy -= p.height;
            if (dy <= -this.halfHeight)
                dy += p.height;

            double dxSquared = Math.pow(dx, 2);
            double dySquared = Math.pow(dy, 2);
            r2 = (float)(dxSquared + dySquared);

            if (r2 < Math.pow(this.particleParticleScreeningLength, 2)) {
                r = (float) Math.sqrt(r2);

                if (r < 0.2) {
//                    System.out.println("WARNING:PARTICLES TOO CLOSE. LOWER CUTOFF FORCE USED");
                    f = 100.0f;
                } else {
                    f = (float) (1 / r2 * Math.exp( -r / this.particleParticleScreeningLength));
                }

//                f = f/r;
                f = f/r + additional_force;

                particles_fx[i] -= f * dx;
                particles_fy[i] -= f * dy;

                particles_fx[j] += f * dx;
                particles_fy[j] += f * dy;
            }


        }
    }

    private void moveParticles() {
        int i;
        float dx, dy;

        for (i = 0; i < nParticles; i++) {
            dx = particles_fx[i] * dt;
            dy = particles_fy[i] * dt;

            particles_x[i] += dx;
            particles_y[i] += dy;

            foldParticleBackPBC(i);

            particles_fx[i] = (float) 0.0;
            particles_fy[i] = (float) 0.0;
        }
    }

    private void foldParticleBackPBC(int i) {
        if (particles_x[i] < 0) {
            particles_x[i] += p.width;
        }

        if (particles_y[i] < 0) {
            particles_y[i] += p.height;
        }

        if (particles_x[i] >= p.width) {
            particles_x[i] -= p.width;
        }

        if (particles_y[i] >= p.height) {
            particles_y[i] -= p.height;
        }
    }

    private void checkVerletRebuildContidionAndSetFlag() {
        int i;
        float dr2;

        this.flagToRebuildVerlet = false;

        for (i = 0; i < this.nParticles; i++) {
            dr2 = (float) (Math.pow(particles_dx_so_far[i], 2) + Math.pow(particles_dy_so_far[i], 2));

            if (dr2 >= verletIntershellSquared) {
                flagToRebuildVerlet = true;
                break;
            }
        }
    }

    private void rebuildVerletList() {
        int i, j;
        float dr2;
        float estimation;

        if (nVerletMax == 0) {
            System.out.println("Verlet list will be built for the first time");

            estimation = nParticles / (float) p.width / (float) p.height;
            System.out.println("System density is " + estimation);

            estimation *= (float) (Math.PI * Math.pow(verletCutoffDistance, 2));
            System.out.println("Particles in R = " + verletCutoffDistance + " shell: " + estimation);

            nVerletMax = (int) estimation * nParticles / 2;
            System.out.println("Estimated nVerletMax: " + nVerletMax);

            this.verletListI = new int[nVerletMax];
            this.verletListJ = new int[nVerletMax];
            System.out.println("OK");
        }

        this.nVerlet = 0;

        for(i = 0; i < this.nParticles; i++) {
            for(j = i + 1; j < this.nParticles; j++) {
                dr2 = distanceSquaredFoldedPBC(particles_x[i], particles_y[i], particles_x[j], particles_y[j]);

                if (dr2 < this.verletCutoffDistanceSquared) {
                    this.verletListI[nVerlet] = i;
                    this.verletListJ[nVerlet] = j;

                    this.nVerlet++;
                    if (nVerlet >= nVerletMax) {
                        System.out.println("Verlet list reallocated from " + nVerletMax);
                        nVerletMax = (int) (1.1 * nVerlet);
                        verletListI = Arrays.copyOf(verletListI, nVerletMax);
                        verletListJ = Arrays.copyOf(verletListI, nVerletMax);
                        System.out.println("New Verlet list max size = " + nVerletMax);
                    }
                }

            }
        }
        System.out.println("Counted Verlet list length: " + nVerlet);
        this.flagToRebuildVerlet = false;

        for (i = 0; i < nParticles; i++) {
            this.particles_dx_so_far[i] = (float) 0.0;
            this.particles_dy_so_far[i] = (float) 0.0;
        }

        for (i = 0; i < nParticles; i++) {
            if (particles_direction[i] == 1) {
                colors[i] = 1;
            } else {
                colors[i] = 0;
            }
        }

        // color particles based on verlet list
//        for (i = 0; i < nVerlet; i++) {
//            if(verletListI[i] == 30) {
//                colors[verletListJ[i]] = 4;
//            }
//
//            if(verletListI[i] == 30) {
//                colors[verletListI[i]] = 4;
//            }
//        }
    }

    private float distanceSquaredFoldedPBC(float x0, float y0, float x1, float y1) {
        float dr2, dx, dy;
        dx = x1 - x0;
        dy = y1 - y0;

        if (dx > this.halfWidth)
            dx -= p.width;
        if (dx <= -this.halfWidth)
            dx += p.width;

        if (dy > this.halfHeight)
            dy -= p.height;
        if (dy <= -this.halfHeight)
            dy += p.height;

        double dxSquared = Math.pow(dx, 2);
        double dySquared = Math.pow(dy, 2);
        dr2 = (float)(dxSquared + dySquared);

        return dr2;
    }

    public void display() {
        for (int i = 0; i < nParticles; i++) {
            if (this.colors[i] == 0) {
                p.stroke(255, 0, 0);
                p.fill(255, 0, 0);
            } else if (this.colors[i] == 1) {
                p.stroke(0, 0, 255);
                p.fill(0, 0, 255);
            } else {
                p.stroke(120, 255, 100);
                p.fill(120, 255, 100);
            }
            p.ellipse(this.particles_x[i], this.particles_y[i], 2, 2);
        }
    }

    public void update() {
        runSimulationStep();
        this.time++;
    }
}
