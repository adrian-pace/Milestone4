package cs211.tangiblegame;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import processing.core.*;
import processing.event.MouseEvent;
import processing.video.*;

public class TangibleGame extends PApplet {
	private final int windowWidth = 640;
	private final int windowHeight = 480;
	private final float boxLength = 300;
	private final float boxHeight = 3;
	private final float ballSize = 13;
	private final float depth = 200;
	private final float gravityConstant = 1.5f;
	private final float cylinderBaseSizeRadius = 18;
	private final float cylinderHeight = 35;
	private final int resolution = 100;
	private final PVector gravity = new PVector(0, gravityConstant, 0);
	private final ArrayList<PVector> cylindersMemory = new ArrayList<PVector>();
	private final Mover mover = new Mover();

	private float speed = 0.5f;
	private boolean shiftPressed = false;
	private float rz = 0;
	private float rx = 0;
	private float ry = 0;
	private float oldMouseX = 0;
	private float oldMouseY = 0;
	private PShape Cylinder;
	private PVector locBeforeCollision;
	private PGraphics dataBackgroundSurface, dataTopViewSurface,
			dataScoreSurface, barChart;
	private float totalScore = 0;
	private float lastScore = 0;
	private ArrayList<Float> scores = new ArrayList<Float>();
	private HScrollbar scrollbar;
	//private PShape tree;
	private PVector rot;
	// private Capture cam;
	private Movie cam;
	private PImage backgroundImage;

	public void setup() {
		size(windowWidth, windowHeight, P3D);
		noStroke();
		// buildCylinder();
		createDataSurfaces();
		Cylinder = loadShape("legotree.obj");

		Cylinder.scale(5);
		Cylinder.rotateX(PI / 2);
		println(Cylinder.height);
		/*
		 * String[] cameras = Capture.list(); if (cameras.length == 0) {
		 * println("There are no cameras available for capture."); exit(); }
		 * else { println("Available cameras:"); for (int i = 0; i <
		 * cameras.length; i++) { println(cameras[i]); } cam = new Capture(this,
		 * cameras[0]); cam.start(); }
		 */
		cam = new Movie(this, "testvideo.mp4"); // Put the video in the same
												// directory

		cam.loop();
		backgroundImage = loadImage("game_background.jpg");
	}

	public void draw() {
		background(backgroundImage);
		// background(200);
		pushMatrix();
		translate(width / 2, height / 2, 0);

		if (!shiftPressed) {
			if (cam.available() == true) {
				cam.read();
			}
			PImage img = cam.get();
			// image(img,-width/2,-height/2);
			pushMatrix();

			if (img.pixels.length > 0) {
				rot = getRotation(img);
				if (rot.y > -PI / 3 && rot.x > -PI / 3 && rot.z > -PI / 3
						&& rot.y < PI / 3 && rot.x < PI / 3 && rot.z < PI / 3) {

					rz = rot.y;
					rx = rot.x;
					ry = rot.z;
				}
			}

			rotateY(ry);
			rotateZ(rz);
			rotateX(rx);
			fill(0, 50, 200, 200);
			box(boxLength, boxHeight, boxLength);
			mover.update();
			mover.checkEdges();
			mover.display();
			drawCylinders();
			mover.checkCollision();

			popMatrix();

		} else {
			fill(0, 50, 200, 200);
			box(boxLength, boxLength, 0);
			for (PVector p : cylindersMemory) {
				pushMatrix();
				translate(p.x, p.y, 0);
				shape(Cylinder);
				popMatrix();
			}
		}
		popMatrix();
		drawAndShowDataSurface();
	}

	void drawCylinders() {
		for (PVector p : cylindersMemory) {
			pushMatrix();
			translate(p.x, 0, p.y);
			rotateX(PI / 2);
			shape(Cylinder);
			popMatrix();
		}
	}

	public void mousePressed() {
		if (!shiftPressed) {
			oldMouseX = mouseX;
			oldMouseY = mouseY;
		} else {
			println(mouseX > -boxLength / 2);
			println(mouseX < boxLength / 2);
			println(mouseY > -boxLength / 2);
			println(mouseY < boxLength / 2);
			float X = mouseX - width / 2;
			float Y = mouseY - height / 2;
			if (X > -boxLength / 2 && X < boxLength / 2 && Y > -boxLength / 2
					&& Y < boxLength / 2) {
				cylindersMemory.add(new PVector(X, Y, 0));
				println(cylindersMemory.size());
			}
		}
	}

	public void mouseDragged() {
		if (!shiftPressed) {
			float diffZ = ((mouseX - oldMouseX) / width) * speed;
			rz = rz + diffZ;
			if (rz > PI / 3) {
				rz = PI / 3;
			} else if (rz < -PI / 3) {
				rz = -PI / 3;
			}
			float diffX = ((mouseY - oldMouseY) / height) * speed;
			rx = rx - diffX;
			if (rx > PI / 3) {
				rx = PI / 3;
			} else if (rx < -PI / 3) {
				rx = -PI / 3;
			}
			oldMouseX = mouseX;
			oldMouseY = mouseY;
		}
	}

	public void mouseWheel(MouseEvent event) {
		float e = event.getCount();
		speed = speed + e * 0.5f;
		if (speed < 0) {
			speed = 0;
		}
		println("THE ACTUAL SPEED IS x" + speed);
	}

	public void keyPressed() {
		if (key == CODED) {
			if (keyCode == LEFT) {
				ry = ry - (PI / 16) * speed;
			} else if (keyCode == RIGHT) {
				ry = ry + (PI / 16) * speed;
			} else if (keyCode == SHIFT) {
				shiftPressed = true;
			}
		}
	}

	public void keyReleased() {
		if (key == CODED) {
			if (keyCode == SHIFT) {
				shiftPressed = false;
			}
		}
	}

	class Mover {
		PVector location;
		PVector velocity;

		Mover() {
			location = new PVector(0, 0, 0);
			velocity = new PVector(0, 0, 0);
		}

		void update() {
			PVector gravityForce = new PVector(sin(rz) * gravityConstant, 0,
					-sin(rx) * gravityConstant);
			float normalForce = 1;
			float mu = 0.01f;
			float frictionMagnitude = normalForce * mu;
			PVector friction = velocity.get();
			friction.mult(-1);
			friction.normalize();
			friction.mult(frictionMagnitude);
			velocity.add(gravityForce);
			velocity.add(friction);
			location.add(velocity);
		}

		void display() {
			pushMatrix();
			translate(location.x, location.y - ballSize, location.z);
			fill(0, 200, 0, 230);
			sphere(ballSize);
			popMatrix();
		}

		void checkEdges() {
			if (location.x > boxLength / 2) {
				velocity.x = velocity.x * -1;
				location.x = boxLength / 2;
				addScore(-velocity.mag());
			}
			if (location.x < -boxLength / 2) {
				velocity.x = velocity.x * -1;
				location.x = -boxLength / 2;
				addScore(-velocity.mag());
			}
			if ((location.z > boxLength / 2)) {
				velocity.z = velocity.z * -1;
				location.z = boxLength / 2;
				addScore(-velocity.mag());
			}
			if (location.z < -boxLength / 2) {
				velocity.z = velocity.z * -1;
				location.z = -boxLength / 2;
				addScore(-velocity.mag());
			}
		}

		void checkCollision() {
			for (PVector p : cylindersMemory) {
				if (sqrt((location.x - p.x) * (location.x - p.x)
						+ (location.z - p.y) * (location.z - p.y)) <= ballSize
						/ 2 + cylinderBaseSizeRadius) {
					PVector locToCyl = new PVector(p.x - location.x, p.y
							- location.z);
					PVector velocity2D = new PVector(velocity.x, velocity.z);
					float teta = PVector.angleBetween(velocity2D, locToCyl);
					if (teta < PI / 2 && teta > -PI / 2) {
						// veloctiy normal to the ball at the collision point
						PVector velocityN = locToCyl.get();
						velocityN.normalize();
						velocityN.mult((velocity.mag() * cos(teta)));
						velocityN.mult(-2);
						PVector velocityN3D = new PVector(velocityN.x, 0,
								velocityN.y);
						velocity.add(velocityN3D);
						location = locBeforeCollision.get();
						addScore(velocity.mag());
						return;
					}
				}
			}
			locBeforeCollision = location.get();
		}
	}

	void buildCylinder() {
		float angle;
		float[] x = new float[resolution + 1];
		float[] y = new float[resolution + 1];
		// get the x and y position on a circle for all the sides
		for (int i = 0; i < x.length; i++) {
			angle = (TWO_PI / resolution) * i;
			x[i] = sin(angle) * cylinderBaseSizeRadius;
			y[i] = cos(angle) * cylinderBaseSizeRadius;
		}
		Cylinder = createShape();
		Cylinder.beginShape(TRIANGLE);
		// draw the border of the cylinder
		for (int i = 0; i < x.length - 1; i++) {
			Cylinder.vertex(x[i], y[i], 0);
			Cylinder.vertex(x[1 + i], y[1 + i], 0);
			Cylinder.vertex(0, 0, 0);
			Cylinder.vertex(x[i], y[i], cylinderHeight);
			Cylinder.vertex(x[1 + i], y[1 + i], cylinderHeight);
			Cylinder.vertex(0, 0, cylinderHeight);
			Cylinder.vertex(x[i], y[i], 0);
			Cylinder.vertex(x[1 + i], y[1 + i], 0);
			Cylinder.vertex(x[i], y[i], cylinderHeight);
			Cylinder.vertex(x[i], y[i], cylinderHeight);
			Cylinder.vertex(x[1 + i], y[1 + i], cylinderHeight);
			Cylinder.vertex(x[1 + i], y[1 + i], 0);
		}
		Cylinder.endShape();
	}

	public void createDataSurfaces() {
		dataBackgroundSurface = createGraphics(width, height / 6, P2D);
		dataTopViewSurface = createGraphics(dataBackgroundSurface.height - 10,
				dataBackgroundSurface.height - 10);
		dataScoreSurface = createGraphics(90, dataBackgroundSurface.height - 10);
		barChart = createGraphics(width - 20 - dataScoreSurface.width
				- dataTopViewSurface.width, dataBackgroundSurface.height - 35);
		scrollbar = new HScrollbar(this, (float) 15 + dataScoreSurface.width
				+ dataTopViewSurface.width, (float) height * 5 / 6 + 10
				+ barChart.height, (float) barChart.width, 20.0f);
	}

	public void drawAndShowDataSurface() {
		dataBackgroundSurface.beginDraw();
		dataBackgroundSurface.background(color(255, 187, 108));
		dataBackgroundSurface.fill(color(255, 187, 108));
		dataBackgroundSurface.noStroke();
		dataBackgroundSurface.rect(0, 0, dataBackgroundSurface.width,
				dataBackgroundSurface.height);
		dataBackgroundSurface.endDraw();
		image(dataBackgroundSurface, 0, height * 5 / 6);

		dataTopViewSurface.noStroke();
		dataTopViewSurface.beginDraw();
		dataTopViewSurface.background(color(51, 51, 255));
		dataTopViewSurface.fill(color(51, 51, 255));
		float ratio = (dataTopViewSurface.width);
		ratio = ratio / boxLength;
		dataTopViewSurface.rect(0, 0, dataTopViewSurface.width,
				dataTopViewSurface.height);
		dataTopViewSurface.fill(255);
		dataTopViewSurface.pushMatrix();
		dataTopViewSurface.translate(dataTopViewSurface.width / 2,
				dataTopViewSurface.height / 2);
		for (PVector p : cylindersMemory) {
			dataTopViewSurface.ellipse(p.x * ratio, p.y * ratio,
					cylinderBaseSizeRadius * ratio, cylinderBaseSizeRadius
							* ratio);
		}
		dataTopViewSurface.fill(0, 200, 0, 230);
		dataTopViewSurface.ellipse(mover.location.x * ratio, mover.location.z
				* ratio, ballSize * ratio, ballSize * ratio);
		dataTopViewSurface.popMatrix();
		dataTopViewSurface.endDraw();
		image(dataTopViewSurface, 5, height * 5 / 6 + 5);

		dataScoreSurface.beginDraw();
		dataScoreSurface.background(color(255, 187, 108));
		dataScoreSurface.noFill();
		dataScoreSurface.stroke(255);
		dataScoreSurface.rect(0, 0, dataScoreSurface.width - 1,
				dataScoreSurface.height - 1);
		float velocity = mover.velocity.mag();
		dataScoreSurface.fill(50);
		dataScoreSurface.text("Total Score:\n" + totalScore + "\n\nVelocity:\n"
				+ velocity + "\n\nLast Score:\n" + lastScore, 5, 5,
				dataScoreSurface.width - 6, dataScoreSurface.height - 6);
		dataScoreSurface.endDraw();
		image(dataScoreSurface, 10 + dataTopViewSurface.width,
				height * 5 / 6 + 5);

		barChart.beginDraw();
		barChart.background(color(255, 187, 108));
		barChart.noFill();
		barChart.stroke(255);
		barChart.rect(0, 0, barChart.width - 1, barChart.height - 1);
		if (scores.size() > 0) {
			int rectWidth = (int) (10 * scrollbar.getPos() + 1);
			int nRectWidth = barChart.width / rectWidth;
			while (scores.size() > nRectWidth) {
				scores.remove(0);
			}
			int rectHeight = 5;
			int nRectHeight = barChart.height / rectHeight;
			int scorePerRect = 100 / nRectHeight;

			barChart.stroke(255);
			fill(color(51, 51, 255));
			for (int i = 0; i < scores.size(); i++) {
				int nRect = (scores.get(i).intValue()) / scorePerRect;
				for (int j = 0; j < nRect; j++) {
					barChart.rect(i * rectWidth, barChart.height - (j + 1)
							* rectHeight, rectWidth, rectHeight);
				}
			}
		}
		barChart.endDraw();
		image(barChart, 15 + dataTopViewSurface.width + dataScoreSurface.width,
				height * 5 / 6 + 5);

		scrollbar.update();
		scrollbar.display();
	}

	public void addScore(float score) {
		lastScore = score;
		totalScore += lastScore;
		if (totalScore < 0) {
			totalScore = 0;
		}
		scores.add(totalScore);
	}

	public PVector getRotation(PImage img) {
		// img = loadImage("board4.jpg");
		// pic 1 115,135
		// pic 2 120, 140
		// pic 3 110, 135
		// pic 4 100, 130
		// room 95,128
		PImage result = hueImage(img, 107, 132);

		result = brightnessThresholding(result, 10, 245);
		// image(result,-width/2,-height/2);
		result = saturationThreshold(result, 100, 255);

		result = convolute(result);
		result = brightnessThresholding(result, 75, 255);

		result = sobel(result, 0.3);

		List<PVector> lines = hough(result, 6);
		QuadGraph quadGraph = new QuadGraph();
		quadGraph.build(lines, img.width, img.height);
		quadGraph.findCycles();
		TwoDThreeD twoDthreeD = new TwoDThreeD(img.width, img.height);
		List<PVector> rotations = new ArrayList<PVector>();
		for (int[] quad : quadGraph.cycles) {
			PVector l1 = lines.get(quad[0]);
			PVector l2 = lines.get(quad[1]);
			PVector l3 = lines.get(quad[2]);
			PVector l4 = lines.get(quad[3]);
			// (intersection() is a simplified version of the
			// intersections() method you wrote last week, that simply
			// return the coordinates of the intersection between 2 lines)

			PVector c12 = intersection(l1, l2);
			PVector c23 = intersection(l2, l3);
			PVector c34 = intersection(l3, l4);
			PVector c41 = intersection(l4, l1);
			if (QuadGraph.isConvex(c12, c23, c34, c41)
					&& QuadGraph.nonFlatQuad(c12, c23, c34, c41)
					&& QuadGraph.validArea(c12, c23, c34, c41, 400000, 10000)) {
				// && QuadGraph.validArea(c12, c23, c34, c41, 2*width*height,
				// 0)) {
				// Choose a random, semi-transparent colour

				List<PVector> list = new ArrayList<PVector>();
				list.add(c12);
				list.add(c23);
				list.add(c34);
				list.add(c41);
				list = CWComparator.sortCorners(list);
				PVector rotationMatrix = twoDthreeD.get3DRotations(list);

				rotations.add(rotationMatrix);

			}
		}
		float rx = 0;
		float ry = 0;
		float rz = 0;
		for (int i = 0; i < rotations.size(); i++) {
			rx = rx + rotations.get(i).x;
			ry = ry + rotations.get(i).y;
			rz = rz + rotations.get(i).z;
		}
		rx = rx / rotations.size();
		ry = ry / rotations.size();
		rz = rz / rotations.size();
		System.out.println("rx=" + rx / (PI * 2) * 360 + ",ry=" + ry / (PI * 2)
				* 360 + ",rz=" + rz / (PI * 2) * 360);
		return new PVector(rx, ry);
	}

	private PVector intersection(PVector line1, PVector line2) {
		double d = Math.cos(line2.y) * Math.sin(line1.y) - Math.cos(line1.y)
				* Math.sin(line2.y);
		int x = (int) Math.round((line2.x * Math.sin(line1.y) - line1.x
				* Math.sin(line2.y))
				/ d);
		int y = (int) Math.round((-line2.x * Math.cos(line1.y) + line1.x
				* Math.cos(line2.y))
				/ d);
		return new PVector(x, y);
	}

	public ArrayList<PVector> hough(PImage edgeImg, int nLines) {
	//	int threshold = 200;
		float discretizationStepsPhi = 0.06f;
		float discretizationStepsR = 2.5f;
		// dimensions of the accumulator
		int phiDim = (int) (Math.PI / discretizationStepsPhi);
		int rDim = (int) (((edgeImg.width + edgeImg.height) * 2 + 1) / discretizationStepsR);
		// our accumulator (with a 1 pix margin around)
		int[] accumulator = new int[(phiDim + 2) * (rDim + 2)];
		// Fill the accumulator: on edge points (ie, white pixels of the edge
		// image), store all possible (r, phi) pairs describing lines going
		// through the point.
		for (int y = 0; y < edgeImg.height; y++) {
			for (int x = 0; x < edgeImg.width; x++) {
				// Are we on an edge?
				if (brightness(edgeImg.pixels[y * edgeImg.width + x]) != 0) {
					// ...determine here all the lines (r, phi) passing through
					// pixel (x,y), convert (r,phi) to coordinates in the
					// accumulator, and increment accordingly the accumulator.
					for (int i = 0; i < phiDim; i++) {
						double phi = i * discretizationStepsPhi;
						double r = x * Math.cos(phi) + y * Math.sin(phi);
						int accPhi = i;
						int accR = (int) (r / discretizationStepsR + (rDim - 1) / 2);
						accumulator[(accPhi + 1) * (rDim + 2) + accR + 1]++;
					}

				}
			}
		}
		PImage houghImg = createImage(rDim + 2, phiDim + 2, ALPHA);
		for (int i = 0; i < accumulator.length; i++) {
			houghImg.pixels[i] = color(min(255, accumulator[i]));
		}
		houghImg.updatePixels();
		// return houghImg;

		ArrayList<Integer> bestCandidates = new ArrayList<Integer>();
		// size of the region we search for a local maximum
		int neighbourhood = 10;
		// only search around lines with more that this amount of votes
		// (to be adapted to your image)
		int minVotes = 150;
		for (int accR = 0; accR < rDim; accR++) {
			for (int accPhi = 0; accPhi < phiDim; accPhi++) {
				// compute current index in the accumulator
				int idx = (accPhi + 1) * (rDim + 2) + accR + 1;
				if (accumulator[idx] > minVotes) {
					boolean bestCandidate = true;
					// iterate over the neighbourhood
					for (int dPhi = -neighbourhood / 2; dPhi < neighbourhood / 2 + 1; dPhi++) {
						// check we are not outside the image
						if (accPhi + dPhi < 0 || accPhi + dPhi >= phiDim)
							continue;
						for (int dR = -neighbourhood / 2; dR < neighbourhood / 2 + 1; dR++) {
							// check we are not outside the image
							if (accR + dR < 0 || accR + dR >= rDim)
								continue;
							int neighbourIdx = (accPhi + dPhi + 1) * (rDim + 2)
									+ accR + dR + 1;
							if (accumulator[idx] < accumulator[neighbourIdx]) {
								// the current idx is not a local maximum!
								bestCandidate = false;
								break;
							}
						}
						if (!bestCandidate)
							break;
					}
					if (bestCandidate) {
						// the current idx *is* a local maximum
						bestCandidates.add(idx);
					}
				}
			}
		}

		Collections.sort(bestCandidates, new HoughComparator(accumulator));
		ArrayList<PVector> returnArray = new ArrayList<PVector>();
		for (int i = 0; i < nLines && i < bestCandidates.size(); i++) {

			// first, compute back the (r, phi) polar coordinates:
			int accPhi = (int) (bestCandidates.get(i) / (rDim + 2)) - 1;
			int accR = bestCandidates.get(i) - (accPhi + 1) * (rDim + 2) - 1;
			float r = (accR - (rDim - 1) * 0.5f) * discretizationStepsR;
			float phi = accPhi * discretizationStepsPhi;
			returnArray.add(new PVector(r, phi));

			// Cartesian equation of a line: y = ax + b
			// in polar, y = (-cos(phi)/sin(phi))x + (r/sin(phi))
			// => y = 0 : x = r / cos(phi)
			// => x = 0 : y = r / sin(phi)
			// compute the intersection of this line with the 4 borders of
			// the image
//			int x0 = 0;
//			int y0 = (int) (r / sin(phi));
//			int x1 = (int) (r / cos(phi));
//			int y1 = 0;
//			int x2 = edgeImg.width;
//			int y2 = (int) (-cos(phi) / sin(phi) * x2 + r / sin(phi));
//			int y3 = edgeImg.width;
//			int x3 = (int) (-(y3 - r / sin(phi)) * (sin(phi) / cos(phi)));

		}
		return returnArray;
	}

	ArrayList<PVector> getIntersections(List<PVector> lines) {
		ArrayList<PVector> intersections = new ArrayList<PVector>();
		for (int i = 0; i < lines.size() - 1; i++) {
			PVector line1 = lines.get(i);
			for (int j = i + 1; j < lines.size(); j++) {
				PVector line2 = lines.get(j);
				// compute the intersection and add it to 'intersections'
				double d = Math.cos(line2.y) * Math.sin(line1.y)
						- Math.cos(line1.y) * Math.sin(line2.y);
				int x = (int) Math.round((line2.x * Math.sin(line1.y) - line1.x
						* Math.sin(line2.y))
						/ d);
				int y = (int) Math
						.round((-line2.x * Math.cos(line1.y) + line1.x
								* Math.cos(line2.y))
								/ d);
				intersections.add(new PVector(x, y));
				// draw the intersection

			}
		}
		return intersections;
	}

	PImage sobel(PImage img, double d) {
		int[][] hKernel = { { 0, 1, 0 }, { 0, 0, 0 }, { 0, -1, 0 } };
		int[][] vKernel = { { 0, 0, 0 }, { 1, 0, -1 }, { 0, 0, 0 } };
		PImage result = createImage(img.width, img.height, ALPHA);
		result.loadPixels();
		img.loadPixels();
		int dim = img.width * img.height;
		float sumh = 0;
		float sumv = 0;
		float[] buffer = new float[dim];
		float max = 0;
		for (int y = 1; y < img.height - 1; y++) {
			for (int x = 1; x < img.width - 1; x++) {
				sumh = 0;
				sumv = 0;
				for (int j = -1; j <= 1; j++) {
					for (int k = -1; k <= 1; k++) {
						sumh += brightness(img.pixels[(y * img.width + x + j
								* img.width + k + dim)
								% dim])
								* hKernel[j + 1][k + 1];
						sumv += brightness(img.pixels[(y * img.width + x + j
								* img.width + k + dim)
								% dim])
								* vKernel[j + 1][k + 1];
						buffer[y * img.width + x] = sqrt(sumh * sumh + sumv
								* sumv);
						if (buffer[y * img.width + x] > max)
							max = buffer[y * img.width + x];
					}

				}

			}

		}
		for (int y = 2; y < img.height - 2; y++) { // Skip top and bottom edges
			for (int x = 2; x < img.width - 2; x++) { // Skip left and right
				if (buffer[y * img.width + x] > (int) (max * d)) {
					result.pixels[y * img.width + x] = color(255);
				} else {
					result.pixels[y * img.width + x] = color(0);
				}
			}
		}
		result.updatePixels();
		return result;
	}

	PImage convolute(PImage img) {
		int[][] matrix = { { 9, 12, 9 }, { 12, 15, 12 }, { 9, 12, 9 } };
		PImage result = createImage(img.width, img.height, ALPHA);
		result.loadPixels();
		img.loadPixels();
		int dim = img.width * img.height;
		int weight = 0;
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				weight += matrix[i][j];
			}
		}
		int sum = 0;
		for (int y = 1; y < img.height - 1; y++) {
			for (int x = 1; x < img.width - 1; x++) {
				for (int j = -1; j <= 1; j++) {
					for (int k = -1; k <= 1; k++) {
						sum += brightness(img.pixels[(y * img.width + x + j
								* img.width + k + dim)
								% dim])
								* matrix[j + 1][k + 1];
					}
				}
				result.pixels[x + y * img.width] = color(sum / weight);
				sum = 0;
			}
		}
		result.updatePixels();
		return result;
	}

	PImage hueImage(PImage img, int lowThreshold, int highThreshold) {
		PImage result = createImage(img.width, img.height, RGB);
		result.loadPixels();
		img.loadPixels();
		for (int i = 0; i < img.pixels.length; i++) {
			int a = img.pixels[i];
			hue(a);
			if (hue(a) <= highThreshold && hue(a) >= lowThreshold)
				result.pixels[i] = img.pixels[i];
			else
				result.pixels[i] = color(0);
		}
		result.updatePixels();
		return result;
	}

	PImage hueImage(PImage img) {
		PImage result = createImage(img.width, img.height, RGB);
		result.loadPixels();
		img.loadPixels();
		for (int i = 0; i < img.width * img.height; i++) {
			result.pixels[i] = color((int) hue(img.pixels[i]));
		}
		result.updatePixels();
		return result;
	}

	PImage brightnessThresholding(PImage img, int lowThreshold,
			int highThreshold) {
		PImage result = createImage(img.width, img.height, RGB);
		result.loadPixels();
		img.loadPixels();
		for (int i = 0; i < img.width * img.height; i++) {
			if (brightness(img.pixels[i]) < lowThreshold
					|| brightness(img.pixels[i]) > highThreshold) {
				result.pixels[i] = color(0);
			} else {
				result.pixels[i] = img.pixels[i];
			}
		}
		result.updatePixels();
		return result;
	}

	PImage saturationThreshold(PImage img, int lowThreshold, int highThreshold) {
		PImage result = createImage(img.width, img.height, RGB);
		result.loadPixels();
		img.loadPixels();
		for (int i = 0; i < img.width * img.height; i++) {
			if (saturation(img.pixels[i]) < lowThreshold
					|| saturation(img.pixels[i]) > highThreshold) {
				result.pixels[i] = color(0);
			} else {
				result.pixels[i] = img.pixels[i];
			}

		}
		result.updatePixels();
		return result;

	}

	PImage intensityThresholding(PImage img, int lowThreshold, int highThreshold) {
		PImage result = createImage(img.width, img.height, RGB);
		result.loadPixels();
		img.loadPixels();
		for (int i = 0; i < img.width * img.height; i++) {
			if ((img.pixels[i]) < lowThreshold
					|| (img.pixels[i]) > highThreshold) {
				result.pixels[i] = color(0);
			} else {
				result.pixels[i] = img.pixels[i];
			}

		}
		result.updatePixels();
		return result;

	}

}
