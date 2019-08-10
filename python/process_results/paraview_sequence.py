from paraview.simple import *

renderView = active_objects.view
source = active_objects.get_source()
displ = GetDisplayProperties()

threshold1 = Threshold(Input=source)

num_frames = source.CellData.GetArray("time").GetRange()[1]

animationScene1 = GetAnimationScene()
animationScene1.NumberOfFrames = int(num_frames)


threshold1ThresholdBetweenTrack = GetAnimationTrack('ThresholdBetween', index=0, proxy=threshold1)

keyFrame9920 = CompositeKeyFrame()

keyFrame9921 = CompositeKeyFrame()
keyFrame9921.KeyTime = 1.0
keyFrame9921.KeyValues = [num_frames]
threshold1ThresholdBetweenTrack.KeyFrames = [keyFrame9920, keyFrame9921]


threshold1ThresholdBetweenTrack_1 = GetAnimationTrack('ThresholdBetween', index=1, proxy=threshold1)

keyFrame8659 = CompositeKeyFrame()
keyFrame8659.KeyTime = 0.0
keyFrame8659.KeyValues = [1.0]

keyFrame8660 = CompositeKeyFrame()
keyFrame8660.KeyTime = 1.0
keyFrame8660.KeyValues = [num_frames + 1]
threshold1ThresholdBetweenTrack_1.KeyFrames = [keyFrame8659, keyFrame8660]

threshold1Display = Show(threshold1, renderView)
Hide(source, renderView)
animationScene1.GoToFirst()
Render()