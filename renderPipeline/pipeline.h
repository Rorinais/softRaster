pipeline::mPipeline = nullptr;
class pipeline {
public:
	pipeline* getInstance() {
		if (mPipeline == nullptr) {
			mPipeline = new pipeline();
		}
		return mPipeline;
	}



private:
	pipeline() = default;
	pipeline(pipeline&) = delete;
	pipeline(pipeline&&) = delete;

	~pipeline() {
		if (mPipeline != nullptr) {
			delete mPipeline;
		}
	}

private:
	static pipeline* mPipeline;

	Raster* mRaster = nullptr;
};
