from itertools import count
import random
import numpy as np


def uniform_dist(elts):
    """
    Uniform distribution over a given finite set of C{elts}
    @param elts: list of any kind of item
    """
    p = 1.0 / len(elts)
    return {e:p for e in elts}


def draw(dist):
    """
    @returns: a randomly drawn element from the distribution
    """
    r = random.random()
    sum = 0.0
    for val, p in dist.items():
        sum += p
        if r < sum:
            return val
    raise Exception('Failed to draw from '+ str(self))


class MDP:
    """
    Needs the following attributes:
    states: list or set of states
    actions: list or set of actions
    discount_factor: real, greater than 0, less than or equal to 1
    start: optional instance of DDist, specifying initial state dist
       if it's unspecified, we will use a uniform over states
    These are functions:
    transition_model: function from (state, action) into distribution over next state
    reward_fn: function from (state, action) to real-valued reward
    """
    
    def __init__(self, states, actions, transition_model, reward_fn, 
                     discount_factor = 1.0, start_dist = None):
        self.states = states
        self.actions = actions
        self.transition_model = transition_model
        self.reward_fn = reward_fn
        self.discount_factor = discount_factor
        self.start = start_dist if start_dist else uniform_dist(states)

    def terminal(self, s):
        """
        Given a state, return True if the state should be considered to
        be terminal.  You can think of a terminal state as generating an
        infinite sequence of zero reward.  Override for your MDP.
        """
        return False

    def init_state(self):
        """
        Randomly choose a state from the initial state distribution
        """
        return draw(self.start)

    def sim_transition(self, s, a):
        """
        Simulate a transition from state s, given action a.  Return
        reward for (s,a) and new state, drawn from transition.  If a
        terminal state is encountered, sample next state from initial
        state distribution
        """
        return (self.reward_fn(s, a),
                self.init_state() if self.terminal(s) else
                    draw(self.transition_model(s, a)))


    def q_update(self, q):
        """
        Perform one step of the value_iteration algorithm and return
        the new q function along with the L1-distance between the
        old and new q-function.
        """
        maxdiff = 0
        q1 = {s:{a:q[s][a] for a in q[s]} for s in q}
        for s in q:
            qs = q[s]
            for a in qs:
                dist = self.transition_model(s, a)
                x = sum(p*max(qs0[a] for a in qs0) for (s0, p) in dist.items() for qs0 in [q[s0]])
                v = self.reward_fn(s, a) + self.discount_factor*x
                maxdiff = max(abs(v - qs[a]), maxdiff)
                q1[s][a] = v
        return q1, maxdiff


    def value_iteration(self, q, eps=0.01, interactive_fn=None):
        """
        Perform value iteration on an MDP.  Terminate when the
        max-norm distance between two successive value function
        estimates is less than eps.
        """
        for h in count():
            q, maxdiff = self.q_update(q)
            if maxdiff < eps:
                break
            if interactive_fn:
                interactive_fn(q, horizon=h)
        return q
    


class TabularQ:
    def __init__(self, states, actions):
        self.actions = actions
        self.states = states
        self.q = dict([((s, a), 0.0) for s in states for a in actions])
    def copy(self):
        q_copy = TabularQ(self.states, self.actions)
        q_copy.q.update(self.q)
        return q_copy
    def set(self, s, a, v):
        self.q[(s,a)] = v
    def get(self, s, a):
        return self.q[(s,a)]
