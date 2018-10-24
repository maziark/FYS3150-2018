TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    vec3.cpp \
    celestialobj.cpp \
    solarsystem.cpp \
    euleralgorithm.cpp \
    verletalgorithm.cpp

HEADERS += \
    vec3.h \
    celestialobj.h \
    solarsystem.h \
    euleralgorithm.h \
    verletalgorithm.h

DISTFILES += \
    initialValues.py
