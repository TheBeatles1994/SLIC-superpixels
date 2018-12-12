QT += core
#QT -= gui
TARGET = SLIC-superpixels
CONFIG += console
CONFIG += app_bundle
CONFIG += c++11

TEMPLATE = app

INCLUDEPATH+=$$PWD/library/win32/opencv/include
             $$PWD/library/win32/opencv/include/opencv
             $$PWD/library/win32/opencv/include/opencv2

SOURCES +=  main.cpp \
    SLIC.cpp \
    glcm.cpp

HEADERS += \
    SLIC.h \
    glcm.h

LIBS+=-L $$PWD/library/win32/opencv/lib/libopencv_*.a
         #$$PWD/library/win32/opencv/bin/libopencv_*.dll


