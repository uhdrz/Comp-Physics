TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    hfcalc.cpp \
    Coulomb_Functions.cpp

HEADERS += \
    hfcalc.h \
    Coulomb_Functions.hpp


LIBS += -llapack -lblas -larmadillo
