#
#    This file is part of WHOI Cable, a program for the static and dynamic
#    analysis of oceanographic cable structures.
#
#    Copyright (C) 2013-2016 by Jason Gobat
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

ifeq ($(OS),windows)
    DIRS = solver model results cli gui
else ifeq ($(OS),nox)
    DIRS = solver model results cli 
else
    DIRS = solver model results widgets cli gui
endif

all::;		@if [ "$(DIRS)" != none ]; then \
		for d in $(DIRS); do \
		echo Making $@ in $(CURDIR)/$$d ...; \
		(cd $$d; make $@) || exit 1; done; fi

clean::;	@if [ "$(DIRS)" != none ]; then \
		for d in $(DIRS); do \
		echo Making $@ in $(CURDIR)/$$d ...; \
		(cd $$d; make $@) || exit 1; done; fi

clobber::;	@if [ "$(DIRS)" != none ]; then \
		for d in $(DIRS); do \
		echo Making $@ in $(CURDIR)/$$d ...; \
		(cd $$d; make $@) || exit 1; done; fi

install::;	@if [ "$(DIRS)" != none ]; then \
		for d in $(DIRS); do \
		if [ -d $$d ]; then \
		echo Making $@ in $(CURDIR)/$$d ...; \
		(cd $$d; make $@) || exit 1; fi; done; fi

depend::;	@if [ "$(DIRS)" != none ]; then \
		for d in $(DIRS); do \
		echo Making $@ in $(CURDIR)/$$d ...; \
		(cd $$d; make $@) || exit 1; done; fi
