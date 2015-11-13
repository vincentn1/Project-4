TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -fopenmp
QMAKE_CXXFLAGS += -O3

LIBS += -fopenmp

INCLUDEPATH += C:\Users\Johannes\Documents\C++\Libraries\

SOURCES += main.cpp \
    lib.cpp \
    random.cpp

HEADERS += \
    lib.h \
    random.h

