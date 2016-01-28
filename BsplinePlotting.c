#include <stdio.h>
#include <stdlib.h>
#include "glut.h"
#include <math.h>



typedef struct {
	double x, y;      // (x, y) coordinates of a 2D point
} Point_2D;


typedef struct {
	int      degree;   // Degree of the B-spline curve
	int      deBoor_count;   // Number of the deBoor points for the B-spline curve
	double   *knots;   // Knot vector of the B-spline curve
	Point_2D  *deBoor_point;     // 2D coordinates of the deBoor points of the B-spline curve
} BSpline;


// Creating a B-spline curve
BSpline  bcr;

// The parameters are used to define a visible window in this application's World Coordinate System. 
double	windowX = 0.0;
double  windowY = 0.0;
double	window_length = 100.0;

int    flag_plotAdaptive = 1;    // A value of 1 indicates Adaptive Rendering and 0 indicates Uniform Rendering
int    flag_dispSampling = 0;    // A value of 1 indicates sampling points will be displayed and 0 indicates thwy will not be displayed
int    flag_dispControlPoly = 1;    // A value of 1 indicates control polygon will be displayed and 0 indicates it will not be displayed
int    sample_count = 10;   // Number of sampling points to be used
double tessEps = 2.0;     // Approximation error for tessellation

static int submenu_id1, submenu_id2, submenu_id3, submenu_id4;
static int menu_id;
static int window;
static int value = 0;



static void Init(void)
{
	glClearColor(1.0, 1.0, 1.0, 0.0);         // To set the color of the display window to white
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(windowX, windowX + window_length, windowY, windowY + window_length);   // To set an orthogonal projection
}

void drawLine(Point_2D a, Point_2D b)
{

	glBegin(GL_LINES);
	glVertex2f(a.x, a.y);
	glVertex2f(b.x, b.y);
	glEnd();

}


//============================================================
void uniformRender()
{
	int h, i, j, k, u, m, deg;
	double t;
	deg = bcr.degree;
	Point_2D* db_point = (Point_2D*)malloc((deg + 1)*sizeof(Point_2D));
	Point_2D* sample_point = (Point_2D*)malloc(sample_count*sizeof(Point_2D));

	for (k = 0; k <= sample_count - 1; k++)
	{
		// Computing sampling intervals within the parameter domain of the curve
		t = (double)(bcr.knots[deg] + (((double)k / (double)(sample_count - 1))*(double)(bcr.knots[bcr.deBoor_count] - bcr.knots[deg])));
		
		j = deg;
		// Determining j such that knots[j]<=t<knots[j+1]
		while (!((t >= bcr.knots[j]) && (t < bcr.knots[j + 1]))&&(j< bcr.deBoor_count))
		{
			j = j + 1;
		}

		//Picking the last segment for the last sampling point
		if (j >= bcr.deBoor_count)
			j = bcr.deBoor_count - 1;


		m = 0;

		//Evaluation of curve at t is impacted by de Boor points [j-degree] to [j]
		for (u = j - deg; u <= j; u++)
		{
			db_point[m].x = bcr.deBoor_point[u].x;
			db_point[m].y = bcr.deBoor_point[u].y;
			m++;

		}
		
		

		// Linear interpolations based on de Boor algorithm
		for (h = 1; h <= deg; h++)
		{
			for (i = j - deg + h, m = 1; i <= j; i++, m++)
			{
				db_point[m - 1].x = (1 - (t - bcr.knots[i]) / (bcr.knots[i + deg + 1 - h] - bcr.knots[i]))*db_point[m - 1].x + ((t - bcr.knots[i]) / (bcr.knots[i + deg + 1 - h] - bcr.knots[i]))*db_point[m].x;
				db_point[m - 1].y = (1 - (t - bcr.knots[i]) / (bcr.knots[i + deg + 1 - h] - bcr.knots[i]))*db_point[m - 1].y + ((t - bcr.knots[i]) / (bcr.knots[i + deg + 1 - h] - bcr.knots[i]))*db_point[m].y;

			}
		}
		sample_point[k] = db_point[0];
	}

	// Joining all sampling points with straight lines
	for (i = 0; i < sample_count - 1; i++)
		drawLine(sample_point[i], sample_point[i + 1]);

	if (flag_dispSampling != 0)
	{
		glPointSize(4.0);            // Display the sampling points
		glBegin(GL_POINTS);
		for (i = 0; i <= sample_count - 1; i++)
			glVertex2f(sample_point[i].x, sample_point[i].y);
		glEnd();
	}


	free(db_point);
	free(sample_point);
}




//============================================================
void extractBezier(Point_2D* bez, int ind)
{
	int     i, j, t;
	int     k;
	double  knots[50]; // Maximum number of knots can be 50
	Point_2D deBoor_point[30]; //Maximum number of deBoor points can be 30

	k = bcr.degree; // k is now the degree of the B-spline curve

					// Curve segment defined over knot span [u(ind),u(ind+1)] is contained in the convex hull of deBoor poitns P(ind-k),...,P(ind)
	for (i = ind - k, j = 0; i <= ind; i++) {
		deBoor_point[j].x = bcr.deBoor_point[i].x;
		deBoor_point[j].y = bcr.deBoor_point[i].y;
		j++;
	}
	// A Bezier curve of degree k will have 2k+2 knots
	for (i = ind - k, j = 0; i <= ind + k + 1; i++) {
		knots[j] = bcr.knots[i];
		j++;
	}

	// Inserting knots to make the left end of the curve become a Bezier end
	while (1)
	{
		for (i = k - 1; i>0; i--)
		{
			if (knots[i] < knots[k]) // Since B-spline knots are always in ascending order, this condition will fail when knots[i]=knots[k]
			{
				j = i;
				break;
			}
			j = 0;
		}

		if (j == 0) break;

		// Updating the control points
		for (i = 0; i<j; i++)
		{
			deBoor_point[i].x = ((knots[k + 1 + i] - knots[k]) / (knots[k + i + 1] - knots[i + 1]))*deBoor_point[i].x
				+ ((knots[k] - knots[i + 1]) / (knots[k + i + 1] - knots[i + 1]))*deBoor_point[i + 1].x;

			deBoor_point[i].y = ((knots[k + 1 + i] - knots[k]) / (knots[k + i + 1] - knots[i + 1]))*deBoor_point[i].y
				+ ((knots[k] - knots[i + 1]) / (knots[k + i + 1] - knots[i + 1]))*deBoor_point[i + 1].y;
		}
		// Updating the knots
		for (i = 0; i<j; i++)
			knots[i] = knots[i + 1];
		knots[j] = knots[k];
	}

	// Inserting knots to make the right end of the curve become a Bezier end
	while (1)
	{
		for (i = k + 2; i< k + k + 1; i++)
		{
			if (knots[i] > knots[k + 1])
			{
				j = i;
				break;
			}
			j = 0;
		}

		if (j == 0) break;

		// Updating the control points
		for (i = k; i >= j - k; i--)
		{
			deBoor_point[i].x = ((knots[k + i] - knots[k + 1]) / (knots[k + i] - knots[i]))*deBoor_point[i - 1].x
				+ ((knots[k + 1] - knots[i]) / (knots[k + i] - knots[i]))*deBoor_point[i].x;
			deBoor_point[i].y = ((knots[k + i] - knots[k + 1]) / (knots[k + i] - knots[i]))*deBoor_point[i - 1].y
				+ ((knots[k + 1] - knots[i]) / (knots[k + i] - knots[i]))*deBoor_point[i].y;
		}
		// Updating the knots
		for (i = k + k + 1; i>j; i--)
			knots[i] = knots[i - 1];
		knots[j] = knots[k + 1];
	}

	// Returning the Bezier control points
	for (i = 0; i< bcr.deBoor_count; i++)
	{
		bez[i].x = deBoor_point[i].x;
		bez[i].y = deBoor_point[i].y;

	}
}

double maxDistance(Point_2D* bez, int deg)
{
	int i;
	double maxheight;
	double *height;
	height = (double*)malloc(30 * sizeof(double));
	// Computing baseline vector and its magnitude
	Point_2D baselineVector;


	baselineVector.x = bez[deg].x - bez[0].x;
	baselineVector.y = bez[deg].y - bez[0].y;
	double baselineMag = (double)sqrt((baselineVector.x*baselineVector.x) + (baselineVector.y*baselineVector.y));

	double crossprod[3];

	double crossprodMag;
	height[0] = 0;

	// Computing height of a intermediate control points from baseline
	Point_2D cpVector; // Vector from first control point to intermediate control point
	for (i = 1; i < deg; i++)
	{
		cpVector.x = bez[i].x - bez[0].x;
		cpVector.y = bez[i].y - bez[0].y;

		// Computing cross product of baseline vector and control point vector
		//z coordinates of baseline and control point vectors are 0
		crossprod[0] = ((baselineVector.y * 0) - (0 * cpVector.y));
		crossprod[1] = -((baselineVector.x * 0) - (0 * cpVector.x));
		crossprod[2] = ((baselineVector.x * cpVector.y) - (baselineVector.y * cpVector.x));
		crossprodMag = (double)sqrt((crossprod[0] * crossprod[0]) + (crossprod[1] * crossprod[1]) + (crossprod[2] * crossprod[2]));
		if (baselineMag != 0)
			height[i] = (double)(crossprodMag / baselineMag);
		else
			height[i] = 0;
	}

	// Finding maximum height of a control point from baseline
	maxheight = 0;
	for (i = 1; i < deg; i++)
	{
		if (maxheight < height[i])
			maxheight = height[i];
	}


	free(height);
	return maxheight;
}


void midSubdivideBezier(Point_2D* bez, int deg, Point_2D* leftBez, Point_2D* rightBez)
{
	Point_2D* deCasteljau_point = (Point_2D*)malloc((deg + 1)*sizeof(Point_2D));
	int i, k, l, r;

	for (i = 0; i <= deg; i++)
	{
		deCasteljau_point[i].x = bez[i].x;
		deCasteljau_point[i].y = bez[i].y;

	}


	l = 0;
	r = deg;


	for (k = 1; k <= deg; k++)
	{
		leftBez[l] = deCasteljau_point[0];
		rightBez[r] = deCasteljau_point[r];

		l = l + 1;
		r = r - 1;
		for (i = 0; i <= deg - k; i++)
		{
			deCasteljau_point[i].x = (0.5* deCasteljau_point[i].x) + (0.5*deCasteljau_point[i + 1].x);
			deCasteljau_point[i].y = (0.5* deCasteljau_point[i].y) + (0.5*deCasteljau_point[i + 1].y);

		}
	}
	leftBez[l] = deCasteljau_point[0];
	rightBez[r] = deCasteljau_point[r];
	free(deCasteljau_point);
}
//Rendering Bezier curves
void plotBezier(Point_2D* bez, int deg)
{
	Point_2D* leftBez = (Point_2D*)malloc((deg + 1)*sizeof(Point_2D));
	Point_2D* rightBez = (Point_2D*)malloc((deg + 1)*sizeof(Point_2D));;
	double height = maxDistance(bez, deg);

	if (height < tessEps)
	{
		drawLine(bez[0], bez[deg]);
		return;
	}
	else
	{
		midSubdivideBezier(bez, deg, leftBez, rightBez);
		plotBezier(leftBez, deg);
		plotBezier(rightBez, deg);
	}
	free(leftBez);
	free(rightBez);
}



//============================================================
void adaptiveRender()
{
	Point_2D*  bez = (Point_2D*)malloc(30 * sizeof(Point_2D)); // assume the degree is not greater than 29.
	int      i;

	for (i = bcr.degree; i< bcr.deBoor_count; i++) // Determining segments between the kth knot and the (n+1)th knot
	{
		if (fabs(bcr.knots[i] - bcr.knots[i + 1]) < 0.00001) // No segment when adjacent knots are equal
			continue;

		extractBezier(bez, i);        // Extract the i-th Bezier curve

		plotBezier(bez, bcr.degree);   // Adaptively plot a Bezier curve 
	}
	free(bez);
}



//============================================================
static void drawCurve(void)
{
	int i;

	glClear(GL_COLOR_BUFFER_BIT);	// clear display window
	glColor3f(1.0, 0.0, 0.0);   // set line segment color to red


								// Draw the control polygon
	glColor3f(1.0, 0.0, 0.0);
	glLineWidth(3.0);
	if (flag_dispControlPoly != 0) {
		glBegin(GL_LINE_STRIP);      // display the control polygon
		for (i = 0; i<bcr.deBoor_count; i++)
			glVertex2f(bcr.deBoor_point[i].x, bcr.deBoor_point[i].y);
		glEnd();

		glPointSize(6.0);            // display the control points
		glBegin(GL_POINTS);
		for (i = 0; i<bcr.deBoor_count; i++)
			glVertex2f(bcr.deBoor_point[i].x, bcr.deBoor_point[i].y);
		glEnd();
	}

	// Draw the curve
	glLineWidth(2.0);
	if (flag_plotAdaptive) {  // plot adaptively
		glColor3f(0.0, 1.0, 0.0);
		adaptiveRender();
	}
	else {  // plot uniformly
		glColor3f(0.0, 0.0, 1.0);
		uniformRender();
	}


	glFlush();		    // process all openGL routines as quickly as possible	         
	glutSwapBuffers();  // swap buffers to display the current frame
}

//============================================================


//============================================================
static void hotkey(unsigned char k, int x, int y)
{
	// Here we are processing keyboard events.
	switch (k)
	{
	case 27:
		free(bcr.deBoor_point);
		free(bcr.knots);
		exit(0);
		break;

		// Toggle plotting the control polygon
	case 'C':
	case 'c':
		flag_dispControlPoly = !flag_dispControlPoly;
		break;

		// Toggle sampling points
	case 'P':
	case 'p':
		flag_dispSampling = !flag_dispSampling;
		break;

		// Toggle adaptive/uniform plotting
	case 'A':
	case 'a':
		flag_plotAdaptive = !flag_plotAdaptive;
		break;

		// Increase tessellation
	case '+':
	case '=':
		if (flag_plotAdaptive)
		{
			tessEps *= 0.7;
			if (tessEps < 0.5)  tessEps = 0.01;
		}
		else
		{
			sample_count += 1;
			if (sample_count > 100) sample_count = 100;
		}
		break;

		// Decrease tessellation
	case '-':
	case '_':
		if (flag_plotAdaptive) {
			tessEps *= 1.4;
			if (tessEps > 50)  tessEps = 100;
		}
		else {
			sample_count -= 1;
			if (sample_count < 2) sample_count = 2;
		}
		break;
	}
}

//============================================================
void chooseWindow()
{
	int    i;
	double left, right, bottom, top;


	left = right = bcr.deBoor_point[0].x;
	for (i = 1; i< bcr.deBoor_count; i++) {
		if (left > bcr.deBoor_point[i].x)  left = bcr.deBoor_point[i].x;
		if (right < bcr.deBoor_point[i].x) right = bcr.deBoor_point[i].x;
	}

	bottom = top = bcr.deBoor_point[0].y;
	for (i = 1; i< bcr.deBoor_count; i++) {
		if (bottom > bcr.deBoor_point[i].y)  bottom = bcr.deBoor_point[i].y;
		if (top < bcr.deBoor_point[i].y) top = bcr.deBoor_point[i].y;
	}

	window_length = top - bottom;
	if (window_length < right - left) window_length = right - left;

	window_length += 100;
	windowY = bottom - 50;
	windowX = left - 50;
}



//============================================================
int readFile(char* filename)
{
	FILE *fp;
	int  i;

	if ((fp = fopen(filename, "r")) == NULL) return 0;  // fail to open the file

	fscanf(fp, "%d%d", &(bcr.degree), &(bcr.deBoor_count));
	bcr.knots = (double *)malloc((bcr.deBoor_count + bcr.degree + 1)*sizeof(double));
	bcr.deBoor_point = (Point_2D *)malloc(bcr.deBoor_count*sizeof(Point_2D));

	for (i = 0; i <= bcr.deBoor_count + bcr.degree; i++)
		fscanf(fp, "%lf", &(bcr.knots[i]));

	for (i = 0; i< bcr.deBoor_count; i++)
		fscanf(fp, "%lf%lf", &(bcr.deBoor_point[i].x), &(bcr.deBoor_point[i].y));
	fclose(fp);

	chooseWindow();

	return 1;
}
void menu(int num) {
	if (num == 0)
	{
		glutDestroyWindow(window);
		free(bcr.deBoor_point);
		free(bcr.knots);
		exit(0);
	}
	else {
		value = num;
		switch (value)
		{

			// Toggle plotting the control polygon
		case 2:
			flag_dispControlPoly = 1;
			break;

		case 3:
			flag_dispControlPoly = 0;
			break;

			// Toggle sampling points
		case 6:
			flag_dispSampling = 1;
			break;

		case 7:
			flag_dispSampling = 0;
			break;

			// Toggle adaptive/uniform plotting
		case 4:
			flag_plotAdaptive = 1;
			break;

		case 5:
			flag_plotAdaptive = 0;
			break;

			// Increase tessellation
		case 8:
			if (flag_plotAdaptive)
			{
				tessEps *= 0.7;
				if (tessEps < 0.5)  tessEps = 0.01;
			}
			else
			{
				sample_count += 1;
				if (sample_count > 100) sample_count = 100;
			}
			break;

			// Decrease tessellation
		case 9:
			if (flag_plotAdaptive) {
				tessEps *= 1.4;
				if (tessEps > 50)  tessEps = 100;
			}
			else {
				sample_count -= 1;
				if (sample_count < 2) sample_count = 2;
			}
			break;
		}
	}
	glutPostRedisplay();
}

void createMenu(void)
{
	submenu_id1 = glutCreateMenu(menu);
	glutAddMenuEntry("Show", 2);
	glutAddMenuEntry("Hide", 3);

	submenu_id2 = glutCreateMenu(menu);
	glutAddMenuEntry("Adaptive", 4);
	glutAddMenuEntry("Uniform", 5);

	submenu_id3 = glutCreateMenu(menu);
	glutAddMenuEntry("Show", 6);
	glutAddMenuEntry("Hide", 7);

	submenu_id4 = glutCreateMenu(menu);
	glutAddMenuEntry("Increase", 8);
	glutAddMenuEntry("Decrease", 9);

	menu_id = glutCreateMenu(menu);
	glutAddSubMenu("Control Polygon", submenu_id1);
	glutAddSubMenu("Rendering Mechanism", submenu_id2);
	glutAddSubMenu("Sampling Points", submenu_id3);
	glutAddSubMenu("Tesselation", submenu_id4);
	glutAddMenuEntry("Quit", 0);

	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

static void idle(void)
{

	drawCurve();
}


//============================================================
void main(int argc, char *argv[])
{
	// load the curve from a file
	char filename[20];

	printf("\n B-spline Curve Plotting \n");
	printf(" ------------------------------------------------------------------------------");
	printf("\n");
	printf("\n Please enter the filename: ");
	scanf("%s", filename);

	if (readFile(filename) == 0) return;

	// help information
	printf("\n");
	printf(" ------------------------------------------------------------------------------");
	printf("\n");
	printf(" ------------------------------------------------------------------------------");
	printf("\n\n\n\n");
	printf(" USER MANUAL:- ");
	printf("\n");
	printf(" ------------------------------------------------------------------------------");
	printf("\n");
	printf(" Adaptive rendering plots curve in Green.\n Uniform rendering plots curve in Blue.\n");
	printf(" ------------------------------------------------------------------------------");
	printf("\n");
	printf(" KEYBOARD:- \n\n");
	printf(" A/a : Toggle adaptive/uniform plotting (Default adaptive)\n");
	printf(" C/c : Toggle plotting the control polygon (Default On)\n");
	printf(" P/p : Toggle sampling points (Default Off)\n");
	printf(" +   : Increase tessellation\n");
	printf(" -   : Decrease tessellation\n");
	printf(" ESC - Quit program\n");
	printf(" ------------------------------------------------------------------------------");
	printf("\n");
	printf(" MOUSE (Right Click):- \n\n");
	printf(" Rendering Mechanism->Adaptive  : Use adaptive plotting (Default)\n");
	printf(" Rendering Mechanism->Uniform  : Use uniform plotting \n");
	printf(" Control Polygon->Show/Hide : Toggle plotting control polygon (Default Show)\n");
	printf(" Sampling Points->Show/Hide : Toggle sampling points (Default Hide)\n");
	printf(" Tesselation->Increase/Decrease   : Increase/decrease tessellation\n");
	printf(" Quit - Quit program\n");
	printf(" ------------------------------------------------------------------------------");
	printf("\n\n");
	printf(" Developed by Aditi Rungta, NTU, November 2015");
	printf("\n");


	// set up graphics window
	glutInit(&argc, argv);                         // initialize GLUT
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);  // set display mode
	glutInitWindowSize(650, 650);                 // set display window width and height
	glutInitWindowPosition(100, 100);             // set top-left display window position
	window = glutCreateWindow("Aditi Rungta - G1502070B");


	Init();                        // execute initialization procedure
	glutIdleFunc(idle);            // enables us to make interaction.
	glutDisplayFunc(drawCurve);    // send graphics to display window

	glutKeyboardFunc(hotkey);
	createMenu();
	glutMainLoop();                // display everything and wait
}
