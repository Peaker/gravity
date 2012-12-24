-- Transliterated from Python from http://users.softlab.ntua.gr/~ttsiod/gravityRK4.html
{-# LANGUAGE TemplateHaskell, NoMonomorphismRestriction #-}

import Control.Applicative
import Control.Lens
import Control.Monad
import Data.Fixed (mod')
import Data.IORef
import Data.List (partition, foldl1')
import Data.Tensor (Vector2(..))
import Graphics.DrawingCombinators ((%%))
import System.Environment (getArgs)
import Vector2
import qualified Control.Monad.Random as Rand
import qualified Data.Foldable as Foldable
import qualified Data.Set as Set
import qualified Graphics.DrawingCombinators as Draw
import qualified Graphics.UI.GLFW as GLFW

windowSize@(Vector2 width height) = Vector2 900 600

-- the density of the planets - used to calculate their mass
-- from their volume (i.e. via their radius)
density = 0.001

-- the gravity coefficient - it's my universe, i can pick whatever i want :-)
gravitystrength = 1e4

-- | Represents position and velocity
data State = State
  { _sPosition :: Vector2 Double
  , _sVelocity :: Vector2 Double
  } deriving (Show, Eq)
makeLenses ''State

-- | Represents velocity and acceleration
data Derivative = Derivative
  { _dVelocity :: Vector2 Double
  , _dAcceleration :: Vector2 Double
  } deriving (Show)
makeLenses ''Derivative

randUnderI upperBound = (`mod` upperBound) <$> Rand.getRandom
randUnder upperBound = (`mod'` upperBound) <$> Rand.getRandom

-- | Class representing a planet. The "_st" member is an instance of
-- "State", carrying the planet's position and velocity - while the
-- "_m" and "_r" members represents the planet's mass and radius
data Planet = Planet
  { _planetRadius :: Double
  , _planetState :: State
  , _planetMass :: Double
  } deriving (Show, Eq)
makeLenses ''Planet

-- | Reversing the setMassFromRadius formula, to calculate radius from
-- mass (used after merging of two planets - mass is added, and new
-- radius is calculated from this)
massFromRadius r = density*4*pi*(r**3)/3
radiusFromMass mass = (3*mass/(density*4*pi))**0.3333
mkPlanetRadius radius state = Planet radius state $ massFromRadius radius
mkPlanetMass mass state = Planet (radiusFromMass mass) state mass

planetPosition = planetState . sPosition
planetVelocity = planetState . sVelocity
planetMomentum p = p^.planetVelocity * pure (p^.planetMass)

-- A nice example of a planet orbiting around our sun :-)
singlePlanet = mkPlanetRadius 1.5 $ State (Vector2 150 300) (Vector2 0 2)
-- otherwise pick a random position and velocity
randomPlanet =
  Rand.evalRandIO $
  mkPlanetRadius 1.5 <$>
  (State
   <$> (fmap fromInteger <$> (Vector2 <$> randUnderI width <*> randUnderI height))
   <*> (subtract 1.5 <$> (Vector2 <$> randUnder 3 <*> randUnder 3)))

-- | Calculate acceleration caused by other planets on this one.
acceleration mass pos p = pure force * distance / pure dr
  where
    distance = p^.planetPosition - pos
    dsq = Foldable.sum $ distance^2
    dr = sqrt dsq
    force | dsq>1e-10 = gravitystrength * mass * p^.planetMass / dsq
          | otherwise = 0

-- | Part of Runge-Kutta method.
mkDerivative planets mass st =
  Derivative (st^.sVelocity) $
  sum $ map (acceleration mass (st^.sPosition)) planets

-- | Part of Runge-Kutta method.
nextDerivative derivative dt planets mass st =
  mkDerivative planets mass $
  st & sPosition +~ derivative^.dVelocity * pure dt
     & sVelocity +~ derivative^.dAcceleration * pure dt

-- | Runge-Kutta 4th order solution to update planet's pos/vel.
updatePlanet planets dt planet = planet
                                 & planetPosition +~ velAvg * pure dt
                                 & planetVelocity +~ accelAvg * pure dt
  where
    mass = planet^.planetMass
    st = planet^.planetState
    a = mkDerivative planets mass st
    b = nextDerivative a (dt * 0.5) planets mass st
    c = nextDerivative b (dt * 0.5) planets mass st
    d = nextDerivative c dt planets mass st
    velAvg = (a^.dVelocity + 2*(b^.dVelocity + c^.dVelocity) + d^.dVelocity) / 6
    accelAvg = (a^.dAcceleration + 2*(b^.dAcceleration + c^.dAcceleration) + d^.dAcceleration) / 6

translateVector2 (Vector2 x y) = Draw.translate (realToFrac x, realToFrac y)
scaleVector2 (Vector2 x y) = Draw.scale (realToFrac x) (realToFrac y)

maxOn f x y
  | f x >= f y = x
  | otherwise = y

positions = go []
  where
    go _ [] = []
    go l (x:xs) = (l, x, xs) : go (x:l) xs

getPlanetCount = do
  args <- map read <$> getArgs
  return $ case args of
    [x] -> x
    _ -> 30

setCallbacks = do
  keysEvents <- newIORef []
  keysPressed <- newIORef Set.empty
  let updateSet True = Set.insert
      updateSet False = Set.delete
      keyCallback key t = do
        modifyIORef keysEvents ((key, t):)
        modifyIORef keysPressed (updateSet t key)
      getEvents = atomicModifyIORef keysEvents $ \events -> ([], reverse events)
  GLFW.setWindowCloseCallback $ fail "Window closed"
  GLFW.setKeyCallback keyCallback
  return (getEvents, keysPressed)

planetsTouch p1 p2 =
  sqrt (Foldable.sum ((p1^.planetPosition - p2^.planetPosition)^2)) <
  p1^.planetRadius + p2^.planetRadius

updatePlanets planetsRef sun dt = modifyIORef planetsRef newPlanets
  where
    -- See if we should merge the ones that are close enough to touch,
    -- using elastic collisions (conservation of total momentum)
    mergeGroup = foldl1' mergeTwo
    mergePlanets [] = []
    mergePlanets (p:ps) = mergeGroup (p : touchers) : mergePlanets others
      where (touchers, others) = partition (planetsTouch p) ps
    mergeTwo p1 p2 =
      mkPlanetMass (p1^.planetMass + p2^.planetMass) State
      { _sPosition = maxOn (^.planetMass) p1 p2 ^. planetPosition
      , _sVelocity =
        (planetMomentum p1 + planetMomentum p2) / pure (p1^.planetMass + p2^.planetMass)
      }
    newPlanets planets =
      mergePlanets . filter (not . planetsTouch sun) $
      -- Calculate the contributions of all the others to its acceleration
      -- (via the gravity force) and update its position and velocity
      map (\(pres, p, posts) -> updatePlanet (sun : pres ++ posts) dt p) $
      positions planets

renderPlanets zoomV shouldClearScreen planets = render image
  where
    planetImage p =
      translateVector2 (p^.planetPosition) %%
      scaleVector2 (pure (p^.planetRadius)) %%
      Draw.circle
    render
      | shouldClearScreen = Draw.clearRender
      | otherwise = Draw.render
    image =
      scaleVector2 zoomV %%
      translateVector2 (-1) %%
      scaleVector2 (2 / fmap fromIntegral windowSize) %%
      Draw.mconcat (map planetImage planets)

handleKeys bClearScreen zoomRef getEvents keysSet = do
  let onKey k = when (Set.member k keysSet)
  -- update zoom factor (numeric keypad +/- keys)
  onKey GLFW.KeyPadAdd      $ modifyIORef zoomRef (* 1.03)
  onKey GLFW.KeyPadSubtract $ modifyIORef zoomRef (/ 1.03)
  onKey GLFW.KeyEsc $ fail "Quit"
  events <- getEvents
  when ((GLFW.CharKey ' ', True) `elem` events) $ do
    modifyIORef bClearScreen not
    newClearScreen <- readIORef bClearScreen
    let verb
          | newClearScreen = "show"
          | otherwise = "hide"
    GLFW.setWindowTitle $
      "Gravity simulation (SPACE: " ++ verb ++ " orbits, keypad +/- : zoom in/out)"

main = do
  planetCount <- getPlanetCount
  True <- GLFW.initialize
  True <- GLFW.openWindow GLFW.defaultDisplayOptions
                          { GLFW.displayOptions_width = width
                          , GLFW.displayOptions_height = height
                          }
  (getEvents, keysSetRef) <- setCallbacks
  -- And God said: Let there be lights in the firmament of the heavens...
  let sun = mkPlanetMass (massFromRadius 1.5 * 1000) State
            { _sPosition = fmap fromIntegral windowSize / 2
            , _sVelocity = 0
            }
  planetsRef <- newIORef . filter (not . planetsTouch sun) =<<
                case planetCount of
                1 -> pure [singlePlanet]
                n -> replicateM n randomPlanet
  -- Zoom factor, changed at runtime via the '+' and '-' numeric keypad keys
  zoomRef <- newIORef (1.0 :: Vector2 Double)
  -- t and dt are unused in this simulation, but are in general,
  -- parameters of engine (acceleration may depend on them)
  let dt = 2
  bClearScreen <- newIORef True
  GLFW.setWindowTitle
    "Gravity simulation (SPACE: show orbits, keypad +/- : zoom in/out)"
  forever $ do
    GLFW.swapBuffers
    shouldClearScreen <- readIORef bClearScreen
    zoomV <- readIORef zoomRef
    renderPlanets zoomV shouldClearScreen . (sun:) =<< readIORef planetsRef
    updatePlanets planetsRef sun dt
    GLFW.pollEvents
    handleKeys bClearScreen zoomRef getEvents =<< readIORef keysSetRef
