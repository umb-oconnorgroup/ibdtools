all: build/ibdtools

build/ibdtools:
	meson build
	ninja -C build ibdtools

build/tests:
	meson build
	ninja -C build tests

ibdtools.AppImage: build/ibdtools

	if [ ! -e ./linuxdeploy-x86_64.AppImage ]; then \
		wget https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage; \
		chmod +x ./linuxdeploy-x86_64.AppImage; \
	fi
	./linuxdeploy-x86_64.AppImage \
		-e build/ibdtools \
		-l $$(ldd build/ibdtools | grep libstdc++ | sed -e 's/^\s*//' | awk '{print $$3}') \
		-l $$(ldd build/ibdtools | grep libhts    | sed -e 's/^\s*//' | awk '{print $$3}') \
		-l $$(ldd build/ibdtools | grep libfmt    | sed -e 's/^\s*//' | awk '{print $$3}') \
		--appdir ibdtools.AppDir \
		-d appimage/ibdtools.desktop \
		-i appimage/icon_256x256.png \
		--output appimage
	rm -rf ibdtools.AppDir
	mv ibdtools*AppImage ibdtools.AppImage

