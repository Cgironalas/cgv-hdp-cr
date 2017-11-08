#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <iostream>

using namespace cv;
using namespace std;

static Mat *masks;
static Mat imageToShow;
static string *filenames;
static int index, folderSize, indexBegin, indexEnd, middle, currentImage;

static void help() {
	cout << "\nThis program demonstrates GrabCut segmentation -- select an object in a region\n"
		"and then grabcut will attempt to segment it out.\n"
		"Call:\n"
		"./grabcut <image_name>\n"
		"\nSelect a rectangular area around the object you want to segment\n" <<
		"\nHot keys: \n"
		"\tESC - quit the program\n"
		"\tr - restore the original image\n"
		"\tn - next iteration\n"
		"\n"
		"\tleft mouse button - set rectangle\n"
		"\n"
		"\tCTRL+left mouse button - set GC_BGD pixels\n"
		"\tSHIFT+left mouse button - set GC_FGD pixels\n"
		"\n"
		"\tCTRL+right mouse button - set GC_PR_BGD pixels\n"
		"\tSHIFT+right mouse button - set GC_PR_FGD pixels\n"
		"\n"
		"\tj - previous image\n"
		"\tk - next image\n"
		"\n"<< endl;
}

const Scalar RED = Scalar(0, 0, 255);
const Scalar PINK = Scalar(230, 130, 255);
const Scalar BLUE = Scalar(255, 0, 0);
const Scalar LIGHTBLUE = Scalar(255, 255, 160);
const Scalar GREEN = Scalar(0, 255, 0);

const int BGD_KEY = EVENT_FLAG_CTRLKEY;
const int FGD_KEY = EVENT_FLAG_SHIFTKEY;

static void getBinMask(const Mat& comMask, Mat& binMask) {
	if (comMask.empty() || comMask.type() != CV_8UC1)
		CV_Error(Error::StsBadArg, "comMask is empty or has incorrect type (not CV_8UC1)");
	if (binMask.empty() || binMask.rows != comMask.rows || binMask.cols != comMask.cols)
		binMask.create(comMask.size(), CV_8UC1);
	binMask = comMask & 1;
}

class GCApplication {
public:
	enum {
		NOT_SET = 0,
		IN_PROCESS = 1,
		SET = 2
	};
	static const int radius = 2;
	static const int thickness = -1;

	bool isInitialized;
	Mat mask;

	void reset();
	void setImageAndWinName(const Mat& _image, const string& _winName);
	void showImage() const;
	void mouseClick(int event, int x, int y, int flags, void* param);
	int nextIter();
	int getIterCount() const { return iterCount; }
	Mat getMask() const { return mask; }

private:
	void setRectInMask();
	void setLblsInMask(int flags, Point p, bool isPr);

	const string* winName;
	const Mat* image;
	Mat bgdModel, fgdModel;
	uchar rectState, lblsState, prLblsState;
	
	Rect rect;
	vector<Point> fgdPxls, bgdPxls, prFgdPxls, prBgdPxls;
	int iterCount;
};


//Set all the elements of the class as null or equivalent (mat with 0s, NOT_SET, 0, etc.)
void GCApplication::reset() {
	if (!mask.empty())
		mask.setTo(Scalar::all(GC_BGD));
	bgdPxls.clear();
	fgdPxls.clear();
	prBgdPxls.clear();
	prFgdPxls.clear();

	isInitialized = false;
	rectState = NOT_SET;
	lblsState = NOT_SET;
	prLblsState = NOT_SET;
	iterCount = 0;
}


//Kind of like a constructor, just sets the image and winName elements of the class
void GCApplication::setImageAndWinName(const Mat& _image, const string& _winName) {
	if (_image.empty() || _winName.empty())
		return;
	image = &_image;
	winName = &_winName;
	mask.create(image->size(), CV_8UC1);
	reset();
}


//
void GCApplication::showImage() const {
	if (image->empty() || winName->empty())
		return;

	Mat res;
	Mat binMask;
	if (!isInitialized)
		image->copyTo(res);
	else {
		getBinMask(mask, binMask);
		image->copyTo(res, binMask);
	}

	vector<Point>::const_iterator it;
	for (it = bgdPxls.begin(); it != bgdPxls.end(); ++it)
		circle(res, *it, radius, BLUE, thickness);
	for (it = fgdPxls.begin(); it != fgdPxls.end(); ++it)
		circle(res, *it, radius, RED, thickness);
	for (it = prBgdPxls.begin(); it != prBgdPxls.end(); ++it)
		circle(res, *it, radius, LIGHTBLUE, thickness);
	for (it = prFgdPxls.begin(); it != prFgdPxls.end(); ++it)
		circle(res, *it, radius, PINK, thickness);

	if (rectState == IN_PROCESS || rectState == SET)
		rectangle(res, Point(rect.x, rect.y), Point(rect.x + rect.width, rect.y + rect.height), GREEN, 2);

	imshow(*winName, res);
}


//
void GCApplication::setRectInMask() {
	CV_Assert(!mask.empty());
	mask.setTo(GC_BGD);
	rect.x = max(0, rect.x);
	rect.y = max(0, rect.y);
	rect.width = min(rect.width, image->cols - rect.x);
	rect.height = min(rect.height, image->rows - rect.y);
	cout << "1" << endl;
	(mask(rect)).setTo(Scalar(GC_PR_FGD));
	cout << "2" << endl;
}


//
void GCApplication::setLblsInMask(int flags, Point p, bool isPr) {
	vector<Point> *bpxls, *fpxls;
	uchar bvalue, fvalue;
	if (!isPr) {
		bpxls = &bgdPxls;
		fpxls = &fgdPxls;
		bvalue = GC_BGD;
		fvalue = GC_FGD;
	} else {
		bpxls = &prBgdPxls;
		fpxls = &prFgdPxls;
		bvalue = GC_PR_BGD;
		fvalue = GC_PR_FGD;
	}
	if (flags & BGD_KEY) {
		bpxls->push_back(p);
		circle(mask, p, radius, bvalue, thickness);
	}
	if (flags & FGD_KEY) {
		fpxls->push_back(p);
		circle(mask, p, radius, fvalue, thickness);
	}
}


//Listener for mouse events
void GCApplication::mouseClick(int event, int x, int y, int flags, void*) {
	// TODO add bad args check
	switch (event) {
		case EVENT_LBUTTONDOWN: {/*set rect or GC_BGD(GC_FGD) labels*/
			bool isb = (flags & BGD_KEY) != 0,
				isf = (flags & FGD_KEY) != 0;
			if (rectState == NOT_SET && !isb && !isf) {
				rectState = IN_PROCESS;
				rect = Rect(x, y, 1, 1);
			}
			if ((isb || isf) && rectState == SET)
				lblsState = IN_PROCESS;
			break;
		}

		case EVENT_RBUTTONDOWN: {/*set GC_PR_BGD(GC_PR_FGD) labels*/
			bool isb = (flags & BGD_KEY) != 0,
				isf = (flags & FGD_KEY) != 0;
			if ((isb || isf) && rectState == SET)
				prLblsState = IN_PROCESS;
			break;
		}

		case EVENT_LBUTTONUP:
			if (rectState == IN_PROCESS) {
				rect = Rect(Point(rect.x, rect.y), Point(x, y));
				rectState = SET;
				cout << "maybe here" << endl;
				setRectInMask();
				cout << "dies here" << endl;
				CV_Assert(bgdPxls.empty() && fgdPxls.empty() && prBgdPxls.empty() && prFgdPxls.empty());
				showImage();
			}
			if (lblsState == IN_PROCESS) {
				setLblsInMask(flags, Point(x, y), false);
				lblsState = SET;
				showImage();
			}
			break;
		case EVENT_RBUTTONUP:
			if (prLblsState == IN_PROCESS) {
				setLblsInMask(flags, Point(x, y), true);
				prLblsState = SET;
				showImage();
			}
			break;

		case EVENT_MOUSEMOVE:
			if (rectState == IN_PROCESS) {
				rect = Rect(Point(rect.x, rect.y), Point(x, y));
				CV_Assert(bgdPxls.empty() && fgdPxls.empty() && prBgdPxls.empty() && prFgdPxls.empty());
				showImage();
			} else if (lblsState == IN_PROCESS) {
				setLblsInMask(flags, Point(x, y), false);
				showImage();
			} else if (prLblsState == IN_PROCESS) {
				setLblsInMask(flags, Point(x, y), true);
				showImage();
			}
			break;
	}
}


//
int GCApplication::nextIter() {
	if (isInitialized)
		grabCut(*image, mask, rect, bgdModel, fgdModel, 1);
	else {
		if (rectState != SET)
			return iterCount;

		if (lblsState == SET || prLblsState == SET)
			grabCut(*image, mask, rect, bgdModel, fgdModel, 1, GC_INIT_WITH_MASK);
		else
			grabCut(*image, mask, rect, bgdModel, fgdModel, 1, GC_INIT_WITH_RECT);
		
		isInitialized = true;
	}
	iterCount++;

	bgdPxls.clear();
	fgdPxls.clear();
	prBgdPxls.clear();
	prFgdPxls.clear();
	return iterCount;
}

GCApplication gcapp;


//Passthrough function for the GCApplication mouseClick
static void on_mouse(int event, int x, int y, int flags, void* param) {
	gcapp.mouseClick(event, x, y, flags, param);
}

Mat showImageWindow() {
	Mat image = imread(filenames[index], 1);// Set to 1 for RGB, with -1 set to RGB with alpha channel
	if (image.empty()) {
		cout << "\n Durn, couldn't read image filename " << filenames[index] << endl;
		return Mat();
	}

	return image;
}

int main(int argc, char** argv) {
	/*	Step 1: Make program load all images from a folder
	Step 2: Run GrabCut, with all user interaction, on a middle frame (if possible the biggest body part), save resulting mask
	Step 3: When user presses "up arrow" or "down arrow" run GrabCut to the next image in the respective possition
	with a copy of the previously generated mask
	Give user option to alter grabcut inputs in that image.
	Save resulting mask. */

	//Input Errors
	if (argc <= 1) {
		cout << "No dataset directory in input." << endl;
		return -1;
	}
	if (argc <= 2) {
		cout << "No image extension in input." << endl;
		return -2;
	}
	if (argc <= 3) {
		cout << "No index begin in input." << endl;
		return -3;
	}
	if (argc <= 4) {
		cout << "No index end in input." << endl;
		return -4;
	}

	//Set input to variables
	string directory, fileExtension;

	cout << "Read variables" << endl;

	directory = argv[1];
	cout << "Directory: " + directory << endl;

	fileExtension = argv[2];
	cout << "File Extension: " + fileExtension << endl;

	sscanf_s(argv[3], "%i", &indexBegin);
	cout << "Index Begin: " + to_string(indexBegin);

	sscanf_s(argv[4], "%i", &indexEnd);
	cout << "Index End: " + to_string(indexEnd);


	help();

	folderSize = indexEnd - indexBegin;
	filenames = new string[folderSize];
	masks = new Mat[folderSize];
	for (int i = 0; i < folderSize; i++) {
		masks[i].setTo(Scalar::all(0));
	}

	for (int i = indexBegin; i < indexEnd; i++) {
		index = i - indexBegin;
		filenames[index] = directory + "/rgb" + std::to_string(i) + "." + fileExtension;
		if (filenames[index].empty()) {
			cout << "\nEmpty filename" << endl;
			return -5;
		}
		//cout << filenames[index] << endl;
	}

	/*
	Test with .rgb
	char * filename = "1001.rgb";
	FILE * f = NULL;
	errno_t err;
	if (err = fopen_s(&f, filename, "rb") != 0) {
	cout << "Error opening file" << endl;
	}
	if (!f) {
	cout << "bad path" << endl;
	}
	char pixels[4096 * 2700];
	fread(pixels, 4096 * 2700, 1, f);
	fclose(f);

	Mat image(4096, 2700, CV_BGR2GRAY, pixels);
	*/

	//string filename = "1893.bmp";
	int middle = (indexEnd - indexBegin) / 2;
	cout << "Middle: ";
	cout << to_string(middle) << endl;

	currentImage = indexBegin + middle;
	index = middle;

	cout << "Current Image: " + to_string(currentImage) << endl;
	cout << "Index: " + to_string(index) + "\n" << endl;


	const string winName = "image";

	namedWindow(winName, WINDOW_AUTOSIZE);
	//Creates a window with name 'image', possible flags:
	//WINDOW_AUTOSIZE (get size of image, cannot be resized)
	//WINDOW_OPENGL (opengl compatible window)

	setMouseCallback(winName, on_mouse, 0);
	//on_mouse is an in between function that they (in openCV in general) to manage mouse interaction
	//They just send the values to their own defined function gcapp.mouseClick

	imageToShow = showImageWindow();

	if (imageToShow.empty()) {
		cout << "image empty" << endl;
		return -6;
	}

	gcapp.setImageAndWinName(imageToShow, winName);
	gcapp.showImage();

	bool go = false;
	for (;;) {
		char c = (char)waitKey(0);
		switch (c) {
		case '\x1b':
			cout << "Exiting ..." << endl;
			goto exit_main;
		case 'r':
			cout << endl;
			gcapp.reset();
			gcapp.showImage();
			break;
		case 'j':
			masks[index] = gcapp.getMask();
			if (currentImage == indexBegin) {
				cout << "Finished this side" << endl;
				currentImage = indexBegin + middle;
				index = middle;
			}
			else {
				currentImage--;
				index--;
				if (masks[index].empty()) {
					go = true;
					masks[index] = gcapp.getMask();
				}
				//masks[index] = gcapp.getMask();
			}

			cout << "Current Image: " + to_string(currentImage) << endl;
			cout << "Index: " + to_string(index) + "\n" << endl;

			//gcapp.reset();
			gcapp.mask = masks[index];
			//gcapp.isInitialized = true;

			imageToShow = showImageWindow();
			gcapp.showImage();
			if (go == true) {
				go = false;
				goto nextIter;
			}
			break;
		case 'k':
			masks[index] = gcapp.getMask();
			if (currentImage == indexEnd) {
				cout << "Finished this side" << endl;
				currentImage = indexBegin + middle;
				index = middle;
			}
			else {
				currentImage++;
				index++;
				if (masks[index].empty()) {
					go = true;
					masks[index] = gcapp.getMask();
				}
				//masks[index] = gcapp.getMask();
			}

			cout << "Current Image: " + to_string(currentImage) << endl;
			cout << "Index: " + to_string(index) + "\n" << endl;

			//gcapp.reset();
			gcapp.mask = masks[index];
			//gcapp.isInitialized = true;

			imageToShow = showImageWindow();
			gcapp.showImage();
			if (go == true) {
				go = false;
				goto nextIter;
			}
			break;
nextIter:
		case 'n':
			int iterCount = gcapp.getIterCount();
			cout << "<" << iterCount << "... ";
			int newIterCount = gcapp.nextIter();
			if (newIterCount > iterCount) {
				gcapp.showImage();
				cout << iterCount << ">" << endl;
			}
			else
				cout << "rect must be determined>" << endl;
			break;
		}
	}
exit_main:
	destroyWindow("image");
	return 0;
}