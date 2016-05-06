#
#    This file is part of WHOI Cable, a program for the static and dynamic
#    analysis of oceanographic cable structures.
#
#    Copyright (C) 1997-2016 by Woods Hole Oceanographic Institution (WHOI)
#    and Jason Gobat
#
#    WHOI Cable is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    WHOI Cable is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with WHOI Cable.  If not, see <http://www.gnu.org/licenses/>.
#

OBJS    = analysis.o environment.o entry.o objects.o combo.o error.o text.o list.o strsep.o control.o plot.o plotutil.o table.o dialog.o duplicates.o pixbuf2mpeg.o
CC      = gcc

ifeq ($(OS),windows)
    GTK_LIBS = -L../../gtk+extra-2.1.2/gtkextra -lgtkextra-win32-2.0 -lz `pkg-config --libs gtk+-2.0` 
    GTK_CFLAGS = `pkg-config --cflags gtk+-2.0` -I../../gtk+extra-2.1.2 
    EXTRA_CFLAGS = -O6 -Wall -funroll-loops -I../include -DZLIB -DWINDOWS
    LDFLAGS = -mwindows
    WINOBJS = resources.o
else ifeq ($(OS),OSX)
    GTK_CFLAGS = -DNATIVE_OSX \
                 -I/Library/Frameworks/Gtk.framework/Headers \
                 -I/Library/Frameworks/GLib.framework/Headers \
                 -I/Library/Frameworks/Cairo.framework/Headers 
    GTK_LIBS      = -lz -framework Gtk -framework GLib -framework Cairo \
                    -lgtkextra-osx-2.0
    EXTRA_CFLAGS = -g -O2 -Wall -funroll-loops -I../include -DZLIB -mmacosx-version-min=10.4
    AVCODEC_FLAGS = 
    AVCODEC_LIBS = 
    LDFLAGS = -g -mmacosx-version-min=10.4
else
    GTK_LIBS = -lz `pkg-config --libs gtk+-2.0 --libs gtkextra-2.0` 
    GTK_CFLAGS = `pkg-config --cflags gtk+-2.0 --cflags gtkextra-2.0` 
    EXTRA_CFLAGS = -g -Wall -funroll-loops -I../include -DZLIB
    AVCODEC_FLAGS = 
    AVCODEC_LIBS = -lavcodec -lavutil
    LDFLAGS = -g
endif

CFLAGS  = $(GTK_CFLAGS) $(EXTRA_CFLAGS) $(AVCODEC_FLAGS)

LIBS    = ../solver/libsolver.a ../model/libmodel.a ../results/libresults.a $(AVCODEC_LIBS) $(GTK_LIBS) -lm

all:    wcable 

wcable: $(OBJS) $(WINOBJS) wcable.o 
		$(CC) -o wcable wcable.o $(OBJS) $(WINOBJS) $(EXTRA_LDFLAGS) $(LDFLAGS) $(LIBS) 
		-rm -f svn_version.h

resources.o: resources.rc bundle/wcable.ico
	windres -i resources.rc -o resources.o 

wcable.o: svn_version.h

svn_version.h:
	/bin/echo -n '#define SVN_Version "' > svn_version.h
	# git describe --tags --long >>svn_version.h
	/bin/echo -n '1' >> svn_version.h
	/bin/echo '"' >> svn_version.h
	/bin/echo -n '#define TIMESTAMP ' >> svn_version.h
	/bin/echo '"'`date +'%a %b %e %H:%M:%S %Z %Y'`'"' >> svn_version.h
	cat svn_version.h

clean:
		rm -f *.o wcable wcable.exe
