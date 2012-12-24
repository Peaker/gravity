module Vector2 where

import Control.Applicative
import Graphics.Rendering.OpenGL.GL (Vector2(..))

-- OpenGL Vector2 missing instances for Num, Fractional, etc :-(
instance Num a => Num (Vector2 a) where
  (+) = liftA2 (+)
  (*) = liftA2 (*)
  negate = fmap negate
  fromInteger = pure . fromInteger
  signum = fmap signum
  abs = fmap abs
instance Fractional a => Fractional (Vector2 a) where
  (/) = liftA2 (/)
  recip = fmap recip
  fromRational = pure . fromRational
