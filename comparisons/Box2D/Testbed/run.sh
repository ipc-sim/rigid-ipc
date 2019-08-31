# A script to run other simulations in different physics simulators.

# Save the directory of this file
THIS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"

## Box2D
BOX2D_DIR=$THIS_DIR/../../../3rdparty/Box2D
if [ ! -d "$BOX2D_DIR" ]; then
    mkdir -p $BOX2D_DIR; git clone git@github.com:erincatto/Box2D.git $BOX2D_DIR
fi
cd $BOX2D_DIR
if [ ! -d "$BOX2D_DIR/Build" ]; then
    premake5 gmake # May need to chage this to `premake gmake`
fi
if [ ! -f "$BOX2D_DIR/Testbed/Tests/TestEntries.cpp.bak" ]; then
    mv $BOX2D_DIR/Testbed/Tests/TestEntries.cpp \
       $BOX2D_DIR/Testbed/Tests/TestEntries.cpp.bak
fi
if [ ! -f "$BOX2D_DIR/Testbed/Framework/Main.cpp.bak" ]; then
    mv $BOX2D_DIR/Testbed/Framework/Main.cpp \
       $BOX2D_DIR/Testbed/Framework/Main.cpp.bak
fi

# Copy tests
for file in $THIS_DIR/Tests/*
do
    file=$(basename -- $file)
    echo "Copying: $file"
    ln -sf $THIS_DIR/Tests/$file $BOX2D_DIR/Testbed/Tests/$file
done

# Copy the main file
ln -sf $THIS_DIR/Framework/Main.cpp $BOX2D_DIR/Testbed/Framework/Main.cpp

make -C Build
if [ $? -eq 0 ]; then
    cd Testbed; ../Build/bin/x86_64/Debug/Testbed
fi
