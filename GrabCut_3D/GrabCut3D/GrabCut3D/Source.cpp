#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <fstream>
#include <iostream>

using namespace cv;
using namespace std;

static bool go;

static Mat imageToShow;

static string sourceDir,	//Directory of the source images
			  maskDir,		//Directory of pre-existing masks or where they will be saved.
			  binMaskDir,	//Directory of where the binary versions of the masks will be saved.
			  fileExtension,
			  tempDir;		//Temporal string for various usages.

static string *imgFiles,	// names of all the images from sourceDir
			  *maskFiles;	// names of all the masks from maskDir

static int index, currentImage, folderSize, indexBegin, middle, indexEnd, firstMask, lastMask;


static void help() {
	cout << "\nThis program demonstrates GrabCut segmentation -- select an object in a region\n"
		"and then grabcut will attempt to segment it out.\n"
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
		"\n"
		"\th - first image with mask\n"
		"\tl - last image with mask\n"
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
	(mask(rect)).setTo(Scalar(GC_PR_FGD));
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
				setRectInMask();
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


void readConfigurationFile(string configFile) {
	ifstream file;
	file.open(configFile);

	if (!file) {
		cerr << "Unable to read configuration file: " + configFile << endl;
		exit(-7);
	}

	int i = 0;
	string curLine;
	while (getline(file, curLine) && i < 6) {
		if (i == 0) {
			sourceDir = curLine.substr(curLine.find(" = ") + 3);
		}
		else if (i == 1) {
			maskDir = curLine.substr(curLine.find(" = ") + 3);
		}
		else if (i == 2) {
			binMaskDir = curLine.substr(curLine.find(" = " ) + 3);
		}
		else if (i == 3) {
			fileExtension = curLine.substr(curLine.find(" = ") + 3);
		}
		else if (i == 4) {
			indexBegin = stoi(curLine.substr(curLine.find(" = ") + 3));
		}
		else if (i == 5) {
			indexEnd = stoi(curLine.substr(curLine.find(" = ") + 3));
		}

		i++;
	}
	if (i < 6) {
		exit(-8);
	}
}

Mat showImageWindow() {
	Mat image = imread(imgFiles[index], 1);// Set to 1 for RGB, with -1 set to RGB with alpha channel
	if (image.empty()) {
		cout << "\n Durn, couldn't read image filename " << imgFiles[index] << endl;
		return Mat();
	}

	return image;
}

string saveCurrentMask() {
	Mat tempMask = gcapp.getMask().clone();

	string maskFile = maskDir + "/rgb" + to_string(currentImage) + ".bmp";
	string binFile = binMaskDir + "/rgb" + to_string(currentImage) + ".bmp";

	imwrite(maskFile, tempMask);

	cout << "Saved mask for image: " + to_string(currentImage) << endl << endl;

	Mat binMask(tempMask.size(), tempMask.type());
	cv::threshold(tempMask, binMask, 2, 255, cv::THRESH_BINARY);

	imwrite(binFile, binMask);


	return maskFile;
}

void resetCurrentImageValues() {
	currentImage = indexBegin + middle;
	index = middle;
}

void checkFollowingMask() {
	if (maskFiles[index] == "") {
		go = true;
		maskFiles[index] = saveCurrentMask();
	}
}

void nextImage() {
	if (currentImage == indexEnd) {
		cout << "Finished right side" << endl;
		resetCurrentImageValues();
	}
	else {
		currentImage++;
		index++;
		checkFollowingMask();
	}
}

void previousImage() {
	if (currentImage == indexBegin) {
		cout << "Finished left side" << endl;
		resetCurrentImageValues();
	}
	else {
		currentImage--;
		index--;
		checkFollowingMask();
	}
}

void followingImage(bool next) {
	maskFiles[index] = saveCurrentMask();

	if (next) {
		nextImage();
	}
	else {
		previousImage();
	}

	cout << "Current Image: " + to_string(currentImage) << endl;
	cout << "Index: " + to_string(index) + "\n" << endl;

	gcapp.mask = imread(maskFiles[index], 0);

	imageToShow = showImageWindow();
	gcapp.showImage();
}

int main(int argc, char** argv) {

	//Input Error
	if (argc <= 1) {
		cout << "No configuration file in input." << endl;
		exit(-1);
	}

	readConfigurationFile(argv[1]);

	cout << endl;
	cout << "Directory:\t" + sourceDir << endl;
	cout << "Mask dir:\t" + maskDir << endl;
	cout << "Bin dir:\t" + binMaskDir << endl;

	cout << "File Extension:\t" + fileExtension << endl;
	cout << "Index Begin:\t" + to_string(indexBegin) << endl;
	cout << "Index End:\t" + to_string(indexEnd) << endl;


	help();


	folderSize = indexEnd - indexBegin + 1;
	cout << "Folder size: " + to_string(folderSize) << endl;
	

	imgFiles = new string[folderSize];
	maskFiles = new string[folderSize];


	for (int i = 0; i < folderSize; i++) {
		maskFiles[i] = "";
	}

	for (int i = 0; i <= folderSize; i++) {
		index = i - indexBegin;
		imgFiles[index] = sourceDir + "/rgb" + std::to_string(i) + "." + fileExtension;
		if (imgFiles[index].empty()) {
			cout << "\nEmpty filename" << endl;
			return -5;
		}
	}

	middle = (indexEnd - indexBegin) / 2;
	cout << "Middle: ";
	cout << to_string(middle) << endl;

	resetCurrentImageValues();

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

		case 'c':

			break;

		case 'l':

			for (int i = folderSize - 1; i >= 0; i--) {

			}
			break;

		case 'h':
			for (int i = 0; i < folderSize; i++) {

			}
			break;

		case 'j':
			followingImage(false);
			if (go == true) {
				go = false;
				goto nextIter;
			}
			break;

		case 'k':
			followingImage(true);
			if (go == true) {
				go = false;
				goto nextIter;
			}
			break;

		case 'f':
			for (int i = 0; i < folderSize; i++) {
				tempDir = maskDir + "rgb" + to_string(i) + ".jpg";
				cout << "Temp Dir: " + tempDir << endl;
				//imwrite(tempDir, masks[i]);
				cout << "Wrote this" << endl;
			}
			goto exit_main;
			break;
nextIter:
		case 'n':
			int iterCount = gcapp.getIterCount();
			cout << "<" << iterCount << "... ";
			int newIterCount = gcapp.nextIter();
			if (newIterCount > iterCount) {
				gcapp.showImage();
				cout << iterCount << ">" << endl;
				maskFiles[index] = saveCurrentMask();
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